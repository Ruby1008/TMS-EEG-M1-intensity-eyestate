function export_preproc_summary(outDir)
% EXPORT_PREPROC_SUMMARY  Collect every Logs\*_procLog.mat into Logs\preproc_summary.json and
% draw the grand-average figure (QC\GrandAverage_ROI.png) for the HTML reports.
data_dir = 'path/to/data';   % release: set this path
if nargin < 1, outDir = [data_dir '\Revision_2026\Preproc']; end
F = dir(fullfile(outDir, 'Logs', '*_procLog.mat'));
S = {};
GA = struct();
for f = F'
    L = load(fullfile(f.folder, f.name)); P = L.procLog;
    d = struct();
    d.id = P.dataset;
    tok = regexp(P.dataset, 'spTMS_(\d+)_(\w+)', 'tokens', 'once');
    d.sub = str2double(tok{1}); d.cond = tok{2};
    d.cnt = P.cntFile; d.started = P.started; d.finished = P.finished;
    d.versions = P.versions; d.params = P.params;
    pi = P.pairInfo;
    d.pair = struct('source', pi.labelSource, 'nPulsesDetected', P.nPulsesDetected, 'nDouble', pi.nDoubleDetections, ...
        'nAfterDedup', pi.nPulses, 'nCodes', pi.nCodes, 'nUnpaired', pi.nUnpaired, 'rho', pi.rho, ...
        'levels', unique(pi.intensity), 'nPer', pi.nPerIntensity, 'medianPulseSize', pi.medianPulseSize);
    if isfield(pi, 'lag_ms'), d.pair.lagMedian = median(pi.lag_ms, 'omitnan'); d.pair.lagP10 = prctile(pi.lag_ms, 10); d.pair.lagP90 = prctile(pi.lag_ms, 90); end
    if isfield(pi, 'triallistRows'), d.pair.triallistRows = pi.triallistRows; end
    d.badChannels = P.badChannels;
    d.trialsEpoched = P.trialsEpoched; d.trialsBefore = P.trialsBefore; d.trialsKept = P.trialsKept;
    d.trialsPerIntensityKept = P.trialsPerIntensityKept;
    d.trialsRejectedPerIntensity = P.trialsRejectedPerIntensity;
    d.ica1 = table2struct(P.ica1); d.ica2 = table2struct(P.ica2);
    d.ica1Reused = P.ica1Reused; d.ica2Reused = P.ica2Reused;
    d.ica1Overrides = P.ica1Overrides; d.ica2Overrides = P.ica2Overrides;
    di = P.decayInfo;
    d.decay = struct('components', di.components, 'groups', {di.groups}, 'ratio', di.ratio, 'peak', di.peak);
    d.blinkElecsPresent = P.blinkElecsPresent; d.moveElecsPresent = P.moveElecsPresent;
    if isfield(di, 'flagA'), d.decay.flagA = di.flagA; d.decay.flagC = di.flagC; end
    d.r2Flagged = struct('ic', {}, 'reason', {}, 'varPct', {}, 'roiFrac', {}, 'topChan', {}, 'topZ', {});
    if isfield(P, 'r2Flagged'), d.r2Flagged = P.r2Flagged; end
    d.pcaDim = []; if isfield(P.params, 'pcaDim'), d.pcaDim = P.params.pcaDim; end
    % trial2023 written into the final file?
    fF = fullfile(outDir, 'Final', [P.dataset '_F_ERM.set']); d.trial2023Filled = NaN;
    if exist(fF, 'file')
        F = load(fF, '-mat'); if isfield(F, 'EEG'), F = F.EEG; end
        if isfield(F.event, 'trial2023'), v = [F.event.trial2023]; d.trial2023Filled = sum(~isnan(v)) / numel(F.event); else, d.trial2023Filled = 0; end
    end
    d.qc = P.tepQC;
    % baseline drift check (README 7): slope of the ROI average over -1000..-100 ms per intensity
    E = P.erp.final; t = E.times; w = t >= -1000 & t <= -100;
    d.baselineDrift_uVper100ms = arrayfun(@(i) 100 * pf1(t(w), E.roi(i, w)), 1:numel(E.levels));
    S{end + 1} = d; %#ok<AGROW>
    key = matlab.lang.makeValidName(d.cond);
    if ~isfield(GA, key), GA.(key) = struct('times', E.times, 'levels', E.levels, 'roi', {{}}); end
    GA.(key).roi{end + 1} = E.roi;
end
fid = fopen(fullfile(outDir, 'Logs', 'preproc_summary.json'), 'w', 'n', 'UTF-8');
fwrite(fid, jsonencode(S, 'PrettyPrint', true, 'ConvertInfAndNaN', true), 'char'); fclose(fid);

% grand averages per condition x intensity (ERM version, ROI)
conds = fieldnames(GA);
if isempty(conds), return, end
fig = figure('Visible', 'off', 'Position', [50 50 560 * numel(conds) 420]);
cols = [0.30 0.55 0.85; 0.95 0.60 0.10; 0.80 0.15 0.15];
for c = 1:numel(conds)
    G = GA.(conds{c}); R = cat(3, G.roi{:});             % level x time x subject
    ax = subplot(1, numel(conds), c); hold(ax, 'on');
    for i = 1:numel(G.levels)
        m = mean(R(i, :, :), 3, 'omitnan');
        plot(ax, G.times, m, 'Color', cols(i, :), 'LineWidth', 1.3, 'DisplayName', sprintf('%d%% RMT', G.levels(i)));
    end
    xlim(ax, [-100 400]); xline(ax, 0, 'k:', 'HandleVisibility', 'off');
    yl = ylim(ax); patch(ax, [-5 10 10 -5], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], 'EdgeColor', 'none', 'HandleVisibility', 'off');
    title(ax, sprintf('%s (n = %d)', strrep(conds{c}, '_', ' '), size(R, 3)));
    xlabel(ax, 'Time from TMS pulse (ms)'); ylabel(ax, 'ROI FC3/FC1/C3/C1 (\muV), 0.1-80 Hz');
    legend(ax, 'Box', 'off', 'Location', 'best'); box(ax, 'off');
end
exportgraphics(fig, fullfile(outDir, 'QC', 'GrandAverage_ROI.png'), 'Resolution', 600); close(fig);
end

function s = pf1(x, y)
ok = ~isnan(y); if sum(ok) < 3, s = NaN; return, end
p = polyfit(x(ok), y(ok), 1); s = p(1);
end
