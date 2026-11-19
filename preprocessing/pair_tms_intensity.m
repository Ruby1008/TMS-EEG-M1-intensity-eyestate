function [EEG, info] = pair_tms_intensity(EEG, evCodes, codeIntensity, varargin)
% PAIR_TMS_INTENSITY  Give every TMS pulse event (from pop_tesa_findpulse) an
% intensity taken from the nearest preceding trigger code, and cross-check it
% against the size of the pulse artefact.
%
% Run on CONTINUOUS data, after pop_tesa_findpulse and before pop_epoch.
%
%   [EEG, info] = pair_tms_intensity(EEG, {'1','2','3'}, [80 100 120], ...
%                   'maxLag', 100, 'minISI', 4000, 'chan', 'FC3', ...
%                   'triallist', triallist, 'triallistCol', 3)
%
% Inputs
%   evCodes        trigger codes, in the same order as codeIntensity
%   codeIntensity  intensity (%RMT) for each code
% Options
%   'tmsLabel'     label of the pulse events            (default 'TMS')
%   'maxLag'       max code-to-pulse delay in ms        (default 100)
%   'minISI'       pulses closer than this (ms) to the previous kept pulse are
%                  treated as double detections and dropped   (default 4000)
%   'chan'         channel used to measure pulse size   (default 'FC3')
%   'triallist'    fallback when the file has no trigger codes at all
%   'triallistCol' column of triallist holding the intensity (default 3)
%
%   'labelSource'  'auto' | 'code' | 'triallist'          (default 'auto')
%   'plotFile'     save the pulse-size boxplot to this file instead of showing it
%
% 2026-10-09 changes (approved by the user): triallist column 3 (1/2/3) is mapped to
% codeIntensity; 'labelSource' forces the triallist for files whose codes are broken
% (e.g. 022/024 Real_EC, where nearly every code is '2').
%
% Output
%   EEG.event(k).intensity  for every kept TMS event (unpaired pulses removed)
%   info                    struct with counts and the pulse-size check

p = inputParser;
p.addParameter('tmsLabel', 'TMS');
p.addParameter('maxLag', 100);
p.addParameter('minISI', 4000);
p.addParameter('chan', 'FC3');
p.addParameter('triallist', []);
p.addParameter('triallistCol', 3);
p.addParameter('labelSource', 'auto');   % 'auto' (codes if present) | 'code' | 'triallist' (2026-10-09: files whose codes are broken)
p.addParameter('plotFile', '');          % save the pulse-size boxplot here (batch mode, figure not shown)
p.addParameter('labelFix', []);          % 2026-10-10: [pulseSize intensity] rows; the pulse whose artefact size matches (<0.5 uV)
                                         % gets this intensity (also keeps pulses that had no code). Source: EMG/VETA cross-check
p.parse(varargin{:});
o = p.Results;

ms2smp = @(ms) ms * EEG.srate / 1000;
types  = cellfun(@num2str, {EEG.event.type}, 'UniformOutput', false);
lats   = [EEG.event.latency];

% --- TMS pulses, double detections removed (keep the first of a burst)
tmsIdx = find(strcmp(types, o.tmsLabel));
[~, ord] = sort(lats(tmsIdx)); tmsIdx = tmsIdx(ord);
keep = true(size(tmsIdx)); lastLat = -Inf;
for k = 1:numel(tmsIdx)
    if lats(tmsIdx(k)) - lastLat < ms2smp(o.minISI)
        keep(k) = false;
    else
        lastLat = lats(tmsIdx(k));
    end
end
info.nDoubleDetections = sum(~keep);
tmsIdx = tmsIdx(keep);

% --- trigger codes
codeIdx = find(ismember(types, evCodes));
info.nPulses = numel(tmsIdx);
info.nCodes  = numel(codeIdx);
intensity = nan(size(tmsIdx));

useTriallist = strcmp(o.labelSource, 'triallist') || (strcmp(o.labelSource, 'auto') && isempty(codeIdx));
if strcmp(o.labelSource, 'triallist') && isempty(o.triallist)
    error('labelSource = triallist but no triallist given.');
end
if ~useTriallist && ~isempty(codeIdx)
    info.labelSource = 'trigger_code';
    codeLat = lats(codeIdx);
    lag = nan(size(tmsIdx));
    for k = 1:numel(tmsIdx)
        d = lats(tmsIdx(k)) - codeLat;              % code must precede the pulse
        d(d < 0 | d > ms2smp(o.maxLag)) = Inf;
        [dmin, j] = min(d);
        if isfinite(dmin)
            intensity(k) = codeIntensity(strcmp(evCodes, types{codeIdx(j)}));
            lag(k) = dmin * 1000 / EEG.srate;
        end
    end
    info.lag_ms = lag;
elseif ~isempty(o.triallist)
    info.labelSource = 'triallist';
    warning('Intensities taken from triallist in pulse order.');
    if size(o.triallist, 1) ~= numel(tmsIdx)
        warning('triallist has %d rows but %d pulses were found - check before use!', ...
            size(o.triallist, 1), numel(tmsIdx));
    end
    n = min(size(o.triallist, 1), numel(tmsIdx));
    tlv = o.triallist(1:n, o.triallistCol);
    % triallist column 3 holds the intensity INDEX 1/2/3 (column 5 is stimulator output), so map it to %RMT
    if all(ismember(tlv, 1:numel(codeIntensity))), tlv = codeIntensity(tlv); end
    intensity(1:n) = tlv;
    info.triallistRows = size(o.triallist, 1);
else
    error('No trigger codes and no triallist given.');
end

% --- manual label fixes (matched on pulse artefact size)
info.labelFixApplied = zeros(0, 3);
if ~isempty(o.labelFix)
    chF = find(strcmpi({EEG.chanlocs.labels}, o.chan));
    ps0 = nan(size(tmsIdx));
    for k = 1:numel(tmsIdx)
        s = round(lats(tmsIdx(k))); w = max(1, s - round(ms2smp(2))) : min(EEG.pnts, s + round(ms2smp(5)));
        ps0(k) = max(abs(diff(double(EEG.data(chF, w)))));
    end
    for r = 1:size(o.labelFix, 1)
        [m, k] = min(abs(ps0 - o.labelFix(r, 1)));
        if m < 0.5
            info.labelFixApplied(end + 1, :) = [k intensity(k) o.labelFix(r, 2)];
            intensity(k) = o.labelFix(r, 2);
        else
            warning('labelFix: no pulse with size %.1f (closest %.1f) - not applied', o.labelFix(r, 1), ps0(k));
        end
    end
end

% --- drop pulses without an intensity
info.nUnpaired = sum(isnan(intensity));
for k = 1:numel(tmsIdx)
    EEG.event(tmsIdx(k)).intensity = intensity(k);
end
dropTms  = setdiff(find(strcmp(types, o.tmsLabel)), tmsIdx(~isnan(intensity)));
EEG.event(dropTms) = [];
EEG = eeg_checkset(EEG, 'eventconsistency');

% --- cross-check: pulse artefact size should grow with intensity
ch = find(strcmpi({EEG.chanlocs.labels}, o.chan));
types = cellfun(@num2str, {EEG.event.type}, 'UniformOutput', false);
tmsIdx = find(strcmp(types, o.tmsLabel));
inten = [EEG.event(tmsIdx).intensity];
psize = nan(size(tmsIdx));
for k = 1:numel(tmsIdx)
    s = round(EEG.event(tmsIdx(k)).latency);
    w = max(1, s - round(ms2smp(2))) : min(EEG.pnts, s + round(ms2smp(5)));
    psize(k) = max(abs(diff(double(EEG.data(ch, w)))));
end
info.intensity = inten;
info.pulseSize = psize;
info.rho = corr(inten(:), psize(:), 'type', 'Spearman');
lev = unique(inten);
info.medianPulseSize = arrayfun(@(v) median(psize(inten == v)), lev);
info.nPerIntensity   = arrayfun(@(v) sum(inten == v), lev);

fprintf('\n[pair_tms_intensity] source: %s | pulses kept: %d | unpaired dropped: %d | double detections: %d\n', ...
    info.labelSource, numel(tmsIdx), info.nUnpaired, info.nDoubleDetections);
for i = 1:numel(lev)
    fprintf('   %3d%% RMT: n = %3d, median pulse size = %.1f\n', lev(i), info.nPerIntensity(i), info.medianPulseSize(i));
end
fprintf('   Spearman rho (pulse size vs intensity) = %.2f\n', info.rho);
if info.rho < 0.3 || any(diff(info.medianPulseSize) <= 0)
    warning('Pulse size does not increase clearly with intensity - check the labels of this dataset!');
end

fig = figure('Name', 'Pulse size vs intensity', 'Visible', fastif_(isempty(o.plotFile), 'on', 'off'));
boxplot(psize, inten); xlabel('Intensity (% RMT)'); ylabel(['Max |diff| at ' o.chan]);
title(sprintf('Spearman \\rho = %.2f (%s)', info.rho, strrep(info.labelSource, '_', ' ')));
if ~isempty(o.plotFile)
    exportgraphics(fig, o.plotFile, 'Resolution', 110); close(fig);
end
end

function v = fastif_(c, a, b)
if c, v = a; else, v = b; end
end
