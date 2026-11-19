function procLog = preprocess_tep(setPrefix, evCodes, triallist, outDir, varargin)
% PREPROCESS_TEP  Pulse-locked TMS-EEG preprocessing for one dataset (EJN revision, 2026).
% Merges TESA_Script_Real_EO_revised.m and TESA_Script_Sham_EO_revised.m into one
% non-interactive function. Parameters follow README section 4 and must not be changed
% without the user's approval; the robustness switches (README 6.7) default to the main analysis.
%
%   procLog = preprocess_tep('spTMS_001_Real_EC', {'1','2','3'}, [], outDir, 'cntFile', cnt)
%
% Required option
%   'cntFile'      raw ANT .cnt file (read only)
% Options (defaults = main analysis)
%   'labelSource'  'auto' | 'code' | 'triallist'      'rate'  findpulse rate (uV/ms), 1000 (user, 2026-10-09)
%   'cutWin' [-5 10]   'baseWin' [-100 -6]   'epochWin' [-2 2]
%   'decayRule' true   'lowpassTEP' 30 ([] = no low-pass)   'causal' false   'reference' 'average' | 'mastoids'
%   'demeanWin' [-1000 -10] (baseline correction right after epoching; user 2026-10-09)
%   'elecNoiseThresh' 6 (ICA round 2; README had 4; user 2026-10-09)
%   'pcaDim' 40, 'decayRuleA' (A u C), 'r2Protect' true  (Step 2, user 2026-10-11; see Review_Step2\Step2_review.html)
%   'overrides'    table from ica_overrides.csv (sub,cond,round,component,action) for this dataset
%   'reuseICA'     true: reuse saved ICA weights when the input data are identical
%   'seed' 1       'saveIntermediate' true
%
% Outputs (outDir): Final\<id>_F_ERM.set / _F_TEP.set, Intermediate\<id>_e ... _eririffiiiRe.set,
% Logs\<id>_pairInfo.mat / _procLog.mat, ICA\<id>_ICA1.png / _ICA2.png / _ICA1_weights.mat / _ICA2_weights.mat /
% _ICA_beforeafter.png, Pairing\<id>_pulsesize.png, QC\<id>_TEP_QC.png / .csv
data_dir = 'path/to/data';   % release: set this path

p = inputParser;
p.addParameter('cntFile', '');
p.addParameter('cedFile', [data_dir '\waveguard64_woRef.ced']);
p.addParameter('codeIntensity', [80 100 120]);
p.addParameter('labelSource', 'auto');
p.addParameter('rate', 1000);
p.addParameter('cutWin', [-5 10]);
p.addParameter('baseWin', [-100 -6]);
p.addParameter('epochWin', [-2 2]);
p.addParameter('decayRule', true);
p.addParameter('lowpassTEP', 30);
p.addParameter('causal', false);
p.addParameter('reference', 'average');
p.addParameter('overrides', []);
p.addParameter('reuseICA', true);
p.addParameter('seed', 1);
p.addParameter('saveIntermediate', true);
p.addParameter('demeanWin', [-1000 -10]);   % added 2026-10-09 (user): baseline-correct each epoch right after epoching
                                            % (raw .cnt carries a ~18 mV DC offset; without this ICA 1 fails). [] = README as written
p.addParameter('elecNoiseThresh', 6);       % changed 2026-10-09 (user): README/TESA default 4 removed components carrying ROI TEP
p.addParameter('labelFix', []);            % [pulseSize intensity] rows from label_overrides.csv (EMG/VETA cross-check, user 2026-10-10)
p.addParameter('decayRuleA', struct('ratioMin', 15, 'latMax', 16, 'halfMax', 40));  % MAIN since 2026-10-11 (user): round-1 decay rule "A u C",
                                            % calibrated on the user's manual review of 6 datasets (Step 2 stage A); [] = original rule C only
p.addParameter('r2Protect', true);          % MAIN since 2026-10-11 (user, Step 2 stage B): keep + flag round-2 comps >10% var whose ROI back-projection >= 30% of the
                                            % (estimated) ROI TEP in any intensity, or whose topography peaks on an ROI electrode (|z| >= 6)
p.addParameter('pcaDim', 40);              % MAIN since 2026-10-11 (user, Step 2 stage C): PCA to 40 dims before ICA round 1 only; [] = no PCA
p.addParameter('trial2023', []);           % 2023/VETA trial number of each paired pulse (time order), from trials_vs_2023_detail.csv
p.addParameter('eFile', '');                % 2026-10-10 (user): start from a 2023 _e dataset (same trials as 2023), epochs re-locked to the pulse, [-1900 1900] ms
p.parse(varargin{:});
o = p.Results;
cutWin = o.cutWin; baseWin = o.baseWin;

dirs = struct('fin', 'Final', 'mid', 'Intermediate', 'log', 'Logs', 'ica', 'ICA', 'pair', 'Pairing', 'qc', 'QC');
for f = fieldnames(dirs)', d = fullfile(outDir, dirs.(f{1})); if ~exist(d, 'dir'), mkdir(d); end, end
D = @(k, name) fullfile(outDir, dirs.(k), name);
saveMid = @(EEG, suffix) saveIfWanted(EEG, [setPrefix suffix], fullfile(outDir, dirs.mid), o.saveIntermediate);

procLog = struct();
procLog.dataset = setPrefix;
procLog.started = char(datetime('now'));
procLog.cntFile = o.cntFile;
procLog.params = rmfield(o, 'overrides');
procLog.params.evCodes = evCodes;
procLog.versions = struct('matlab', version, 'eeglab', eeg_getversion, 'tesa', tesaVersion());
procLog.overridesApplied = o.overrides;

if isempty(o.eFile)
%% Load raw data and channel locations
EEG = pop_loadeep_v4(o.cntFile);
EEG.setname = setPrefix;
EEG = pop_chanedit(EEG, 'load', {o.cedFile, 'filetype', 'autodetect'});
EEG = eeg_checkset(EEG);
ced = readtable(o.cedFile, 'FileType', 'text', 'Delimiter', '\t');
if ~isequal(strtrim(string(ced.labels))', string({EEG.chanlocs.labels}))
    error('%s: channel labels in the .cnt differ from %s - stop and check.', setPrefix, o.cedFile);
end

%% TMS pulses -> intensity
EEG = pop_tesa_findpulse(EEG, 'FC3', 'refract', 3, 'rate', o.rate, 'tmsLabel', 'TMS', 'plots', 'off');
procLog.nPulsesDetected = sum(strcmp(cellfun(@num2str, {EEG.event.type}, 'uni', 0), 'TMS'));
[EEG, pairInfo] = pair_tms_intensity(EEG, evCodes, o.codeIntensity, 'chan', 'FC3', 'maxLag', 100, 'minISI', 4000, ...
    'triallist', triallist, 'labelSource', o.labelSource, 'plotFile', D('pair', [setPrefix '_pulsesize.png']), 'labelFix', o.labelFix);
save(D('log', [setPrefix '_pairInfo.mat']), 'pairInfo');
procLog.pairInfo = pairInfo;
% 2023 / VETA trial number on every TMS event (for single-trial TEP-MEP analyses; user 2026-10-10)
tmsK = find(strcmp(cellfun(@num2str, {EEG.event.type}, 'uni', 0), 'TMS'));
[~, ordK] = sort([EEG.event(tmsK).latency]); tmsK = tmsK(ordK);
t23 = nan(1, numel(tmsK));
if numel(o.trial2023) == numel(tmsK), t23 = o.trial2023(:)';
elseif ~isempty(o.trial2023), warning('%s: trial2023 has %d entries but %d pulses - set to NaN', setPrefix, numel(o.trial2023), numel(tmsK)); end
for k = 1:numel(tmsK), EEG.event(tmsK(k)).trial2023 = t23(k); end
procLog.trial2023 = t23;

%% Epoch on the pulse, remove the pulse window
EEG = pop_epoch(EEG, {'TMS'}, o.epochWin, 'newname', [setPrefix '_e'], 'epochinfo', 'yes');
EEG = eeg_checkset(EEG);
procLog.trialsEpoched = EEG.trials;
if ~isempty(o.demeanWin)
    EEG = pop_rmbase(EEG, o.demeanWin);   % removes the ~18 mV amplifier DC offset before the zero-filled cut and ICA 1
end
else
%% Start from the 2023 _e epochs (read only), re-locked to the pulse
[EEG, pairInfo] = epochs_from_e(o.eFile, evCodes, o.codeIntensity, 'triallist', triallist, ...
    'labelSource', fastifStr(strcmp(o.labelSource, 'triallist'), 'triallist', 'code'), 'rate', o.rate, 'winMs', [-1900 1900]);
ced = readtable(o.cedFile, 'FileType', 'text', 'Delimiter', '\t');
if ~isequal(strtrim(string(ced.labels))', string({EEG.chanlocs.labels}))
    error('%s: channel labels in the _e file differ from %s - stop and check.', setPrefix, o.cedFile);
end
save(D('log', [setPrefix '_pairInfo.mat']), 'pairInfo');
procLog.pairInfo = pairInfo;
procLog.nPulsesDetected = pairInfo.nPulses;
procLog.trialsEpoched = EEG.trials;
fig = figure('Visible', 'off'); boxplot(pairInfo.pulseSize, pairInfo.intensity);
xlabel('Intensity (% RMT)'); ylabel('Max |diff| at FC3'); title(sprintf('Spearman \\rho = %.2f (%s)', pairInfo.rho, strrep(pairInfo.labelSource, '_', ' ')));
exportgraphics(fig, D('pair', [setPrefix '_pulsesize.png']), 'Resolution', 110); close(fig);
if ~isempty(o.demeanWin)
    EEG = pop_rmbase(EEG, o.demeanWin);
end
end
saveMid(EEG, '_e');
EEG = pop_tesa_removedata(EEG, cutWin, [], {'TMS'});
saveMid(EEG, '_er');

%% Bad channels (baseline log-variance, robust z > 5) and bad trials (joint probability 5 SD)
isEOG       = contains(upper({EEG.chanlocs.labels}), 'EOG');
eegChanlocs = EEG.chanlocs(~isEOG);
bl  = EEG.times >= -1000 & EEG.times <= -10;
lv  = median(log(squeeze(var(EEG.data(:, bl, :), 0, 2))), 2);
z   = (lv - median(lv(~isEOG))) / (1.4826 * mad(lv(~isEOG), 1));
badChan = {EEG.chanlocs(abs(z) > 5 & ~isEOG').labels};
procLog.badChanZ = table(string({EEG.chanlocs.labels})', z, 'VariableNames', {'channel', 'z'});
if ~isempty(badChan), EEG = pop_select(EEG, 'nochannel', badChan); end
nTrialsBefore = EEG.trials;
intenBefore = epochIntensity(EEG);
EEG = pop_jointprob(EEG, 1, 1:EEG.nbchan, 5, 5, 0, 1, 0, [], 0);
EEG = eeg_checkset(EEG);
procLog.badChannels  = badChan;
procLog.trialsBefore = nTrialsBefore;
procLog.trialsAfterJointprob = EEG.trials;
procLog.trialsRejectedPerIntensity = countPer(intenBefore, o.codeIntensity) - countPer(epochIntensity(EEG), o.codeIntensity);
fprintf('[%s] bad channels: %s | trials kept: %d of %d\n', setPrefix, strjoin(badChan, ' '), EEG.trials, nTrialsBefore);
erp.preICA1 = summaryERP(EEG, o.codeIntensity);

%% ICA round 1: decay rule + TESA tmsMuscle (classification on the full component set, removal once)
[EEG, w1] = runICA(EEG, D('ica', [setPrefix '_ICA1_weights.mat']), o.seed, o.reuseICA, o.pcaDim);
procLog.ica1Reused = w1.reused;
if o.decayRule
    [~, decayInfo] = decay_component_rule(EEG, 'win', [13 30], 'preWin', [-500 -10], 'ratio', 25, 'minPeak', 10, 'ruleA', o.decayRuleA, ...
        'remove', false, 'plot', false);
else
    decayInfo = struct('components', [], 'ratio', [], 'peak', [], 'groups', {{}}, 'params', 'decay rule off');
end
EEG = pop_tesa_compselect(EEG, 'compCheck', 'off', 'remove', 'off', 'saveWeights', 'off', ...
    'tmsMuscle', 'on', 'tmsMuscleThresh', 8, 'tmsMuscleWin', [11 30], 'tmsMuscleFeedback', 'off', ...
    'blink', 'off', 'move', 'off', 'muscle', 'off', 'elecNoise', 'off');
cls1 = lastTesaClass(EEG);
reasons = repmat({''}, 1, size(EEG.icaweights, 1));
reasons(decayInfo.components) = {'decay'};
for c = cls1.rejectTmsMuscle, reasons{c} = joinReason(reasons{c}, 'tmsMuscle'); end
[reasons, ovr1] = applyOverrides(reasons, o.overrides, 1);
ica1 = componentTable(EEG, reasons, cutWin);
ica1.tmsMuscleRatio = tmsMuscleRatio(EEG, [11 30]);
procLog.decayInfo = decayInfo;
procLog.ica1 = ica1;
procLog.ica1Overrides = ovr1;
procLog.icaCompClass1 = cls1;
icaOverviewFigure(EEG, ica1, cutWin, sprintf('%s  ICA round 1 (removed %d of %d)', setPrefix, sum(ica1.removed), height(ica1)), ...
    D('ica', [setPrefix '_ICA1.png']), decayInfo);
rem1 = find(ica1.removed)';
if ~isempty(rem1), EEG = pop_subcomp(EEG, rem1, 0); end
EEG = eeg_checkset(EEG);
erp.postICA1 = summaryERP(EEG, o.codeIntensity);
saveMid(EEG, '_eri');

%% Interpolate, filter
EEG = pop_tesa_removedata(EEG, cutWin, [], {'TMS'});
saveMid(EEG, '_erir');
EEG = pop_tesa_interpdata(EEG, 'cubic', [20 20]);
saveMid(EEG, '_eriri');
EEG = butterFilter(EEG, 0.1, 80, 'bandpass', o.causal);       % 0.1-80 Hz, 4th order
saveMid(EEG, '_eririf');
EEG = butterFilter(EEG, 59, 61, 'bandstop', o.causal);        % 59-61 Hz, 4th order
saveMid(EEG, '_eririff');
EEG = pop_tesa_interpdata(EEG, 'cubic', [20 20]);
saveMid(EEG, '_eririffi');
erp.preICA2 = summaryERP(EEG, o.codeIntensity);

%% ICA round 2: blink, eye movement, persistent muscle, electrode noise (tmsMuscle OFF)
[EEG, w2] = runICA(EEG, D('ica', [setPrefix '_ICA2_weights.mat']), o.seed, o.reuseICA, []);
procLog.ica2Reused = w2.reused;
present = {EEG.chanlocs.labels};
% blink / eye-movement electrodes: README Fp1/Fp2 and F7/F8; fall back to the nearest available pair only when
% both are missing (bad channels), e.g. spTMS_023_Sham_EO (2026-10-12)
blinkE = intersect({'Fp1', 'Fp2'}, present, 'stable');
if isempty(blinkE), blinkE = intersect({'Fpz', 'AF3', 'AF4'}, present, 'stable'); end
moveE = intersect({'F7', 'F8'}, present, 'stable');
if isempty(moveE), moveE = intersect({'F5', 'F6'}, present, 'stable'); end
if isempty(moveE), moveE = intersect({'FT7', 'FT8'}, present, 'stable'); end
procLog.blinkElecsPresent = blinkE;
procLog.moveElecsPresent  = moveE;
EEG = pop_tesa_compselect(EEG, 'compCheck', 'off', 'remove', 'off', 'saveWeights', 'off', ...
    'tmsMuscle', 'off', ...
    'blink', 'on', 'blinkThresh', 2.5, 'blinkElecs', blinkE, 'blinkFeedback', 'off', ...
    'move', 'on', 'moveThresh', 2, 'moveElecs', moveE, 'moveFeedback', 'off', ...
    'muscle', 'on', 'muscleThresh', -0.31, 'muscleFreqIn', [7 70], 'muscleFreqEx', [48 52], 'muscleFeedback', 'off', ...
    'elecNoise', 'on', 'elecNoiseThresh', o.elecNoiseThresh, 'elecNoiseFeedback', 'off');
cls2 = lastTesaClass(EEG);
reasons = repmat({''}, 1, size(EEG.icaweights, 1));
for c = cls2.rejectBlink,          reasons{c} = joinReason(reasons{c}, 'blink'); end
for c = cls2.rejectEyeMove,        reasons{c} = joinReason(reasons{c}, 'eyeMove'); end
for c = cls2.rejectMuscle,         reasons{c} = joinReason(reasons{c}, 'muscle'); end
for c = cls2.rejectElectrodeNoise, reasons{c} = joinReason(reasons{c}, 'elecNoise'); end
r2flag = struct('ic', {}, 'reason', {}, 'varPct', {}, 'roiFrac', {}, 'topChan', {}, 'topZ', {});
if o.r2Protect
    [reasons, r2flag] = protectR2(EEG, reasons, o.codeIntensity, baseWin);
end
procLog.r2Flagged = r2flag;
[reasons, ovr2] = applyOverrides(reasons, o.overrides, 2);
ica2 = componentTable(EEG, reasons, cutWin);
procLog.ica2 = ica2;
procLog.ica2Overrides = ovr2;
procLog.icaCompClass2 = cls2;
icaOverviewFigure(EEG, ica2, cutWin, sprintf('%s  ICA round 2 (removed %d of %d)', setPrefix, sum(ica2.removed), height(ica2)), ...
    D('ica', [setPrefix '_ICA2.png']), []);
rem2 = find(ica2.removed)';
if ~isempty(rem2), EEG = pop_subcomp(EEG, rem2, 0); end
EEG = eeg_checkset(EEG);
erp.postICA2 = summaryERP(EEG, o.codeIntensity);
saveMid(EEG, '_eririffii');

%% Linear interpolation of the pulse window, drop EOG, interpolate bad channels, re-reference, baseline
EEG = pop_tesa_interpdata(EEG, 'linear');
saveMid(EEG, '_eririffiii');
eogIdx = find(contains(upper({EEG.chanlocs.labels}), 'EOG'));
if ~isempty(eogIdx), EEG = pop_select(EEG, 'nochannel', eogIdx); EEG = eeg_checkset(EEG); end
if ~isempty(badChan), EEG = pop_interp(EEG, eegChanlocs, 'spherical'); EEG = eeg_checkset(EEG); end
switch o.reference
    case 'average'
        EEG = pop_reref(EEG, []);
    case 'mastoids'
        EEG = pop_reref(EEG, find(ismember({EEG.chanlocs.labels}, {'M1', 'M2'})), 'keepref', 'on');
    otherwise
        error('unknown reference %s', o.reference);
end
EEG = eeg_checkset(EEG);
saveMid(EEG, '_eririffiiiRe');
EEG = pop_rmbase(EEG, baseWin);

%% Output 1: ERM (0.1-80 Hz + 59-61 Hz band-stop); QC on this version (R1-2)
EEG_ERM = EEG;
tepQC = tep_subject_qc(EEG_ERM, setPrefix, outDir, cutWin);
erp.final = summaryERP(EEG_ERM, o.codeIntensity);
beforeAfterFigure(erp, cutWin, setPrefix, D('ica', [setPrefix '_ICA_beforeafter.png']));

%% Output 2: TEP (additional 30 Hz low-pass)
if isempty(o.lowpassTEP)
    EEG_TEP = EEG;
else
    EEG_TEP = butterFilter(EEG, [], o.lowpassTEP, 'lowpass', o.causal);
end

%% Event labels '80'/'100'/'120' and save
for v = 1:2
    if v == 1, E = EEG_ERM; suffix = '_F_ERM'; else, E = EEG_TEP; suffix = '_F_TEP'; end
    E.event = E.event(strcmp(cellfun(@num2str, {E.event.type}, 'uni', 0), 'TMS'));
    for i = 1:numel(E.event), E.event(i).type = num2str(E.event(i).intensity); end
    E = eeg_checkset(E, 'eventconsistency');
    E.setname = [setPrefix suffix];
    pop_saveset(E, 'filename', [setPrefix suffix '.set'], 'filepath', fullfile(outDir, dirs.fin));
end

%% Log
procLog.trialsKept = EEG.trials;
procLog.trialsPerIntensityPaired = pairInfo.nPerIntensity;
procLog.trialsPerIntensityKept   = countPer(epochIntensity(EEG), o.codeIntensity);
procLog.tepQC = tepQC;
procLog.erp = erp;
procLog.finished = char(datetime('now'));
save(D('log', [setPrefix '_procLog.mat']), 'procLog');
end

%% ---------------------------------------------------------------- helpers
function saveIfWanted(EEG, name, folder, doSave)
if doSave
    EEG.setname = name;
    pop_saveset(EEG, 'filename', [name '.set'], 'filepath', folder);
end
end

function v = tesaVersion()
v = 'unknown';
f = which('eegplugin_tesa');
if ~isempty(f), tok = regexp(f, 'TESA([\d.]+)', 'tokens', 'once'); if ~isempty(tok), v = tok{1}; end, end
end

function inten = epochIntensity(EEG)
% intensity of the time-0 TMS event of each epoch
inten = nan(1, EEG.trials);
for k = 1:numel(EEG.event)
    ev = EEG.event(k);
    if strcmp(num2str(ev.type), 'TMS') && isfield(ev, 'intensity') && ~isempty(ev.intensity)
        lat = mod(ev.latency - 1, EEG.pnts) + 1;
        if abs(EEG.times(round(lat))) < 1, inten(ev.epoch) = ev.intensity; end
    end
end
end

function n = countPer(inten, levels)
n = arrayfun(@(v) sum(inten == v), levels);
end

function r = joinReason(r, new)
if isempty(r), r = new; else, r = [r '+' new]; end
end

function cls = lastTesaClass(EEG)
f = fieldnames(EEG.icaCompClass);
f = f(strncmp(f, 'TESA', 4));
cls = EEG.icaCompClass.(f{end});
end

function [EEG, info] = runICA(EEG, wfile, seed, reuse, pcaDim)
% FastICA (symm, tanh) with a fixed seed; weights saved so a re-run (e.g. after an override in
% another round) gives identical components. Weights are reused only if the input data are identical.
sig = dataSignature(EEG);
info.reused = false;
if reuse && exist(wfile, 'file')
    W = load(wfile);
    if isequal(W.signature, sig)
        EEG.icaweights = W.icaweights; EEG.icasphere = W.icasphere;
        EEG.icawinv = W.icawinv; EEG.icachansind = W.icachansind; EEG.icaact = [];
        EEG = eeg_checkset(EEG, 'ica');
        info.reused = true;
        fprintf('ICA weights reused from %s\n', wfile);
        return
    end
end
rng(seed);
if isempty(pcaDim)
    EEG = pop_tesa_fastica(EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off');
else
    % same FastICA settings as pop_tesa_fastica, with PCA reduction to pcaDim dimensions, then TESA's sorting
    EEG = pop_runica(EEG, 'icatype', 'fastica', 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off', 'lastEig', pcaDim);
    EEG = tesa_sortcomps(EEG);
    fprintf('FastICA after PCA to %d dimensions.\n', pcaDim);
end
icaweights = EEG.icaweights; icasphere = EEG.icasphere; icawinv = EEG.icawinv; icachansind = EEG.icachansind; %#ok<NASGU>
chanlabels = {EEG.chanlocs.labels}; signature = sig; rngSeed = seed; %#ok<NASGU>
save(wfile, 'icaweights', 'icasphere', 'icawinv', 'icachansind', 'chanlabels', 'signature', 'rngSeed');
end

function s = dataSignature(EEG)
d = double(EEG.data);
s = [size(d, 1) size(d, 2) size(d, 3) sum(d(:)) sum(abs(d(:))) sum(d(:, 1:7:end).^2, 'all')];
s = round(s, 6, 'significant');
end

function [reasons, applied] = applyOverrides(reasons, ovr, round_)
applied = {};
if isempty(ovr), return, end
ovr = ovr(ovr.round == round_, :);
for i = 1:height(ovr)
    c = ovr.component(i);
    if c < 1 || c > numel(reasons), warning('override IC%d outside 1..%d - ignored', c, numel(reasons)); continue, end
    switch lower(strtrim(char(ovr.action(i))))
        case 'remove'   % strip any keep wrapper (flaggedKeep(...) from the round-2 rule, manualKeep(...)) so the component is removed
            r0 = regexprep(reasons{c}, '^(flaggedKeep|manualKeep)\((.*)\)$', '$2');
            reasons{c} = joinReason(r0, 'manual');
        case 'keep',   reasons{c} = ['manualKeep(' reasons{c} ')'];
        otherwise, warning('unknown override action %s', ovr.action(i)); continue
    end
    applied{end + 1} = sprintf('R%d IC%d %s', round_, c, ovr.action(i)); %#ok<AGROW>
end
end

function T = componentTable(EEG, reasons, cutWin)
% one row per component: reason, removed flag, explained variance of the back-projected data
act = eeg_getdatact(EEG, 'component', 1:size(EEG.icaweights, 1));
n = size(act, 1);
e = zeros(n, 1);
for i = 1:n, e(i) = sum(EEG.icawinv(:, i).^2) * mean(act(i, :).^2); end
varPct = 100 * e / sum(e);
avgVar = arrayfun(@(i) var(mean(act(i, :, :), 3)), (1:n)');      % as TESA (variance of the trial average)
tesaVarPct = 100 * avgVar / sum(avgVar);
reasons = reasons(:);
removed = ~cellfun(@isempty, reasons) & ~startsWith(reasons, 'manualKeep') & ~startsWith(reasons, 'flaggedKeep');
T = table((1:n)', string(reasons), removed, varPct, tesaVarPct, ...
    'VariableNames', {'component', 'reason', 'removed', 'varPct', 'tesaVarPct'});
T.Properties.UserData = struct('cutWin', cutWin);
end

function [reasons, fl] = protectR2(EEG, reasons, levels, baseWin)
% Round-2 protection rule (Step 2): among components TESA classified as elecNoise or muscle (not blink/eyeMove; user 2026-10-11)
% that explain > 10% variance, keep (and flag)
% those whose ROI back-projection (average reference, EOG excluded) has RMS 15-200 ms >= 30% of the estimated ROI TEP
% in any intensity, or whose topography peaks on an ROI electrode with max |z| >= 6. The estimated TEP is the ROI average
% after removing every component TESA classified, average-referenced, baseline-corrected (an approximation of the final TEP).
roi = {'FC3', 'FC1', 'C3', 'C1'};
fl = struct('ic', {}, 'reason', {}, 'varPct', {}, 'roiFrac', {}, 'topChan', {}, 'topZ', {});
act = eeg_getdatact(EEG, 'component', 1:size(EEG.icaweights, 1));
n = size(act, 1); e = zeros(n, 1);
for i = 1:n, e(i) = sum(EEG.icawinv(:, i).^2) * mean(act(i, :).^2); end
varPct = 100 * e / sum(e);
lab = {EEG.chanlocs(EEG.icachansind).labels}; nonEOG = ~contains(upper(lab), 'EOG'); ri = ismember(upper(lab), upper(roi));
t = EEG.times; bw = t >= baseWin(1) & t <= baseWin(2); tw = t >= 15 & t <= 200;
inten = epochIntensity(EEG);
rem = find(~cellfun(@isempty, reasons));
X = double(EEG.data(EEG.icachansind, :, :));
clean = X - reshape(EEG.icawinv(:, rem) * reshape(act(rem, :, :), numel(rem), []), size(X));
for c = rem
    if varPct(c) <= 10, continue, end
    % 2026-10-11 (user, option A): only elecNoise / muscle components are protected; blink / eyeMove are always removed
    if ~(contains(reasons{c}, 'elecNoise') || contains(reasons{c}, 'muscle')) || contains(reasons{c}, 'blink') || contains(reasons{c}, 'eyeMove'), continue, end
    frac = zeros(1, numel(levels));
    for i = 1:numel(levels)
        tr = inten == levels(i);
        y = mean(clean(:, :, tr), 3); y = y - mean(y(nonEOG, :), 1); y = mean(y(ri, :), 1); y = y - mean(y(bw));
        x = EEG.icawinv(:, c) * mean(act(c, :, tr), 3); x = x - mean(x(nonEOG, :), 1); x = mean(x(ri, :), 1); x = x - mean(x(bw));
        frac(i) = rms(x(tw)) / max(rms(y(tw)), eps);
    end
    z = zscore(EEG.icawinv(:, c)); [zm, k] = max(abs(z)); ch = lab{k};
    if max(frac) >= 0.3 || (zm >= 6 && any(strcmpi(ch, roi)))
        fl(end + 1) = struct('ic', c, 'reason', reasons{c}, 'varPct', varPct(c), 'roiFrac', frac, 'topChan', ch, 'topZ', zm); %#ok<AGROW>
        reasons{c} = ['flaggedKeep(' reasons{c} ')'];
    end
end
if ~isempty(fl), fprintf('[protectR2] kept and flagged: %s\n', mat2str([fl.ic])); end
end

function r = tmsMuscleRatio(EEG, win)
% TESA tmsMuscle score: mean |trial-average| in the window / mean over the whole epoch
act = eeg_getdatact(EEG, 'component', 1:size(EEG.icaweights, 1));
a = abs(mean(act, 3));
w = EEG.times >= win(1) & EEG.times <= win(2);
r = mean(a(:, w), 2) ./ mean(a, 2);
end

function EEG = butterFilter(EEG, lo, hi, type, causal)
% Zero-phase: TESA pop_tesa_filtbutter (4th order). Causal (robustness only): same design, filter().
if ~causal
    EEG = pop_tesa_filtbutter(EEG, lo, hi, 4, type);
    return
end
fs = EEG.srate;
switch type
    case 'bandpass', [z, p, k] = butter(2, [lo hi] / (fs / 2), 'bandpass');
    case 'bandstop', [z, p, k] = butter(2, [lo hi] / (fs / 2), 'stop');
    case 'lowpass',  [z, p, k] = butter(4, hi / (fs / 2), 'low');
end
sos = zp2sos(z, p, k);
d = double(EEG.data);
for c = 1:size(d, 1)
    x = reshape(d(c, :, :), 1, []);
    d(c, :, :) = reshape(sosfilt(sos, x), 1, EEG.pnts, []);
end
EEG.data = d;
end

function S = summaryERP(EEG, levels)
% ROI (FC3/FC1/C3/C1) average per intensity and EEG-channel butterfly (all trials), for the before/after figure
roiIx = find(ismember(upper({EEG.chanlocs.labels}), upper({'FC3', 'FC1', 'C3', 'C1'})));
eegIx = find(~contains(upper({EEG.chanlocs.labels}), 'EOG'));
inten = epochIntensity(EEG);
S.times = EEG.times;
S.levels = levels;
S.roi = nan(numel(levels), EEG.pnts);
for i = 1:numel(levels)
    tr = inten == levels(i);
    if any(tr), S.roi(i, :) = mean(mean(EEG.data(roiIx, :, tr), 3), 1); end
end
S.butterfly = single(mean(EEG.data(eegIx, :, :), 3));
S.labels = {EEG.chanlocs(eegIx).labels};
end

function beforeAfterFigure(erp, cutWin, setPrefix, file)
st = {'preICA1', 'postICA1', 'preICA2', 'postICA2'};
ttl = {'Before ICA 1 (raw, cut window zeroed)', 'After ICA 1', 'Before ICA 2 (filtered, interpolated)', 'After ICA 2'};
fig = figure('Visible', 'off', 'Position', [50 50 1900 760]);
cols = [0.30 0.55 0.85; 0.95 0.60 0.10; 0.80 0.15 0.15];
for s = 1:4
    S = erp.(st{s}); t = S.times; win = t >= -200 & t <= 500; show = win & ~(t >= cutWin(1) & t <= cutWin(2));
    ax = subplot(2, 4, s); hold(ax, 'on');
    for i = 1:numel(S.levels), plot(ax, t, S.roi(i, :), 'Color', cols(i, :), 'LineWidth', 1.1, 'DisplayName', sprintf('%d%%', S.levels(i))); end
    setLimits(ax, S.roi(:, show), cutWin); title(ax, ttl{s}, 'FontSize', 9);
    if s == 1, ylabel(ax, 'ROI FC3/FC1/C3/C1 (\muV)'); legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7); end
    ax = subplot(2, 4, 4 + s); hold(ax, 'on');
    plot(ax, t, S.butterfly', 'Color', [0.55 0.55 0.55 0.6], 'LineWidth', 0.5);
    setLimits(ax, S.butterfly(:, show), cutWin); xlabel(ax, 'Time from TMS pulse (ms)');
    if s == 1, ylabel(ax, 'Butterfly, all trials (\muV)'); end
end
sgtitle(strrep(setPrefix, '_', '\_'));
exportgraphics(fig, file, 'Resolution', 100); close(fig);
end

function setLimits(ax, y, cutWin)
y = y(isfinite(y)); if isempty(y), y = [-1 1]; end   % e.g. all ROI channels bad at this stage
lim = prctile(y(:), [0.5 99.5]); r = diff(lim); if ~(r > 0), r = 1; end
ylim(ax, lim + [-0.1 0.1] * r); xlim(ax, [-200 500]);
patch(ax, [cutWin fliplr(cutWin)], [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)], [0.85 0.85 0.85], 'EdgeColor', 'none', 'HandleVisibility', 'off');
xline(ax, 0, 'k:', 'HandleVisibility', 'off'); box(ax, 'off');
end

function v = fastifStr(c, a, b)
if c, v = a; else, v = b; end
end
