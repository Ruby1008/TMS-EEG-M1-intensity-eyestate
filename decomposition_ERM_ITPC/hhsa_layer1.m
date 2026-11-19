function info = hhsa_layer1(id, varargin)
% HHSA_LAYER1  .set -> SPM -> first-layer improved CEEMDAN (eemd1_tfpar) for one dataset.
%   hhsa_layer1('spTMS_001_Real_EO')          writes HHSA\spm_<id>.mat and timf_/timph_/timfm_/timam_ files
%   hhsa_layer1(id, 'force', true)            re-run even if timf_ exists
% Same steps as the original (spm_eeg_convert -> spm_eeg_prep defaulteegsens/coregister -> eemd1_tfpar)
% with the parameters in hhsa_setup, condition labels = TMS intensity (80/100/120) from EEG.event,
% mastoids dropped, and rng(P.seed) fixed before the decomposition.
p = inputParser; p.addParameter('force', false); p.parse(varargin{:});
P = hhsa_setup(); spm('defaults', 'eeg'); eeglab nogui; close all
info = struct('id', id, 'nTrials', NaN, 'nChan', NaN, 'minutes', NaN, 'skipped', false);
spmF = fullfile(P.out, ['spm_' id '.mat']); timfF = fullfile(P.out, ['timf_spm_' id '.mat']);
if exist(timfF, 'file') && ~p.Results.force, info.skipped = true; fprintf('%s: timf exists, skip\n', id); return, end
t0 = tic;
% --- 1. EEGLAB -> SPM (61 scalp channels, labels = intensity)
EEG = pop_loadset('filename', [id P.suffix '.set'], 'filepath', P.src);
EEG = pop_select(EEG, 'nochannel', P.dropChan);
inten = nan(1, EEG.trials);
for k = 1:numel(EEG.event)
    if ~isempty(EEG.event(k).intensity), inten(EEG.event(k).epoch) = EEG.event(k).intensity; end
end
assert(~any(isnan(inten)), '%s: epochs without intensity', id);
tmpSet = fullfile(P.out, ['tmp_' id '.set']);
pop_saveset(EEG, 'filename', ['tmp_' id '.set'], 'filepath', P.out, 'savemode', 'onefile');
S = struct('dataset', tmpSet, 'outfile', spmF, 'channels', 'all', 'timewin', [], 'blocksize', 3276800, ...
    'checkboundary', 1, 'eventpadding', 0, 'saveorigheader', 0, 'inputformat', [], 'mode', 'epoched');
S.conditionlabels = arrayfun(@(v) sprintf('%d', v), inten, 'UniformOutput', false);
D = spm_eeg_convert(S);
delete(tmpSet); delete(strrep(tmpSet, '.set', '.fdt'));
D = chantype(D, 1:D.nchannels, 'EEG');
D = conditions(D, 1:D.ntrials, S.conditionlabels); save(D);          % spm_eeg_convert kept only the first label for all epochs
assert(isequal(sort(unique(D.conditions)), sort(unique(S.conditionlabels))), '%s: condition labels not set', id);
S = struct('D', D, 'task', 'defaulteegsens', 'save', 1); D = spm_eeg_prep(S);
% 'coregister' (head-mesh display) is skipped: it opens a 3-D figure (rotate3d) and fails in -batch;
% it only adds D.inv for source analysis, which the sensor-level CEEMDAN/ITPC/HHCFPC steps do not use.
info.nTrials = D.ntrials; info.nChan = D.nchannels;
fprintf('%s: %d channels, %d trials (%s)\n', id, D.nchannels, D.ntrials, strjoin(unique(D.conditions), '/'));
% --- 2. first-layer CEEMDAN
if P.l1.parallel && isempty(gcp('nocreate')), parpool('local', P.l1.CoreN); end
S = P.l1; S.D = D;
rng(P.seed);                                   % noise bank is drawn once, before the trial/channel loops
eemd1_tfpar(S);
info.minutes = toc(t0) / 60;
fprintf('%s: layer 1 done in %.1f min\n', id, info.minutes);
end
