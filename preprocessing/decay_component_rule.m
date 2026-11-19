function [EEG, decayInfo] = decay_component_rule(EEG, varargin)
% DECAY_COMPONENT_RULE  Flag and remove TMS decay-artefact ICA components that the
% TESA tmsMuscle rule misses (TESA compares the early window with the WHOLE-EPOCH
% mean, so slow decays and DC shifts are kept).
%
% Run in ICA round 1, right after pop_tesa_fastica and BEFORE pop_tesa_compselect.
%
% Rule (baseline-referenced, per component, on the trial-averaged |activation|):
%   mean |act| in 13-30 ms  >=  ratio x  mean |act| in the pre-stimulus window
%   AND peak |act| in 13-30 ms >= minPeak
% evaluated on all trials and, separately, within each intensity
% (a component is flagged if it passes in all trials OR in any single intensity).
%
% Options
%   'win'      early window (ms)                    default [13 30]
%   'preWin'   pre-stimulus reference window (ms)   default [-500 -10]
%   'ratio'    threshold ratio                      default 25
%   'minPeak'  minimum peak (activation units)      default 10
%   'remove'   remove flagged components            default true
%   'plot'     show flagged components for review   default true

p = inputParser;
p.addParameter('win', [13 30]);
p.addParameter('preWin', [-500 -10]);
p.addParameter('ratio', 25);
p.addParameter('minPeak', 10);
p.addParameter('remove', true);
p.addParameter('plot', true);
p.addParameter('ruleA', []);   % 2026-10-10 (user, Step 2 stage A): struct('ratioMin',15,'latMax',16,'halfMax',40). A component is ALSO
                               % flagged if ratio >= ratioMin in any group (no peak threshold), the all-trial |average| peaks at
                               % <= latMax ms (searched 10-60 ms) and falls below half of its 11-ms value by halfMax ms ("A u C")
p.parse(varargin{:});
o = p.Results;

if isempty(EEG.icaact)
    EEG.icaact = eeg_getdatact(EEG, 'component', 1:size(EEG.icaweights, 1));
end
act = EEG.icaact;                                  % comp x time x trial
t   = EEG.times;
w   = t >= o.win(1)    & t <= o.win(2);
pw  = t >= o.preWin(1) & t <= o.preWin(2);

% intensity of each epoch (from pair_tms_intensity)
inten = nan(1, EEG.trials);
if isfield(EEG.event, 'intensity')
    for k = 1:numel(EEG.event)
        if ~isempty(EEG.event(k).intensity), inten(EEG.event(k).epoch) = EEG.event(k).intensity; end
    end
end
groups = {true(1, EEG.trials)};
lev = unique(inten(~isnan(inten)));
for v = lev, groups{end + 1} = inten == v; end %#ok<AGROW>

nComp   = size(act, 1);
flag    = false(nComp, 1);
ratioAll = zeros(nComp, numel(groups));
peakAll  = zeros(nComp, numel(groups));
for g = 1:numel(groups)
    a = abs(mean(act(:, :, groups{g}), 3));        % evoked activation
    ratioAll(:, g) = mean(a(:, w), 2) ./ max(mean(a(:, pw), 2), eps);
    peakAll(:, g)  = max(a(:, w), [], 2);
    flag = flag | (ratioAll(:, g) >= o.ratio & peakAll(:, g) >= o.minPeak);
end

flagC = flag;
decayInfo.flagC = flagC';
if ~isempty(o.ruleA)
    a = abs(mean(act, 3)); lw = find(t >= 10 & t <= 60);
    [~, k] = max(a(:, lw), [], 2); peakLat = t(lw(k))';
    i0 = find(t >= 11, 1); halfT = nan(nComp, 1);
    for c = 1:nComp
        j = find(a(c, i0:end) < a(c, i0) / 2, 1);
        if ~isempty(j), halfT(c) = t(i0 + j - 1); end
    end
    flagA = any(ratioAll >= o.ruleA.ratioMin, 2) & peakLat <= o.ruleA.latMax & halfT <= o.ruleA.halfMax;
    flag = flag | flagA;
    decayInfo.flagA = flagA'; decayInfo.peakLatency_ms = peakLat'; decayInfo.halfDecay_ms = halfT';
    fprintf('[decay_component_rule] rule A flags %s; original rule (C) flags %s\n', mat2str(find(flagA)'), mat2str(find(flagC)'));
end
decayInfo.components = find(flag)';
decayInfo.ratio      = ratioAll;
decayInfo.peak       = peakAll;
decayInfo.groups     = [{'all'}, arrayfun(@(v) sprintf('%d%%', v), lev, 'UniformOutput', false)];
decayInfo.params     = o;
fprintf('[decay_component_rule] flagged %d component(s): %s\n', numel(decayInfo.components), mat2str(decayInfo.components));

if o.plot && any(flag)
    idx = find(flag);
    figure('Name', 'Decay components (round 1)');
    for i = 1:numel(idx)
        subplot(numel(idx), 2, 2 * i - 1);
        plot(t, mean(act(idx(i), :, :), 3), 'k'); xlim([-200 500]); xline(0, ':');
        title(sprintf('IC%d  ratio %.0f  peak %.0f', idx(i), max(ratioAll(idx(i), :)), max(peakAll(idx(i), :))));
        subplot(numel(idx), 2, 2 * i);
        topoplot(EEG.icawinv(:, idx(i)), EEG.chanlocs(EEG.icachansind), 'electrodes', 'off');
    end
end

if o.remove && any(flag)
    EEG = pop_subcomp(EEG, find(flag), 0);
    EEG = eeg_checkset(EEG);
end
end
