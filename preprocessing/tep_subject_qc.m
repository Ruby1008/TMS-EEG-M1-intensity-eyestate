function qc = tep_subject_qc(EEG, setPrefix, outDir, cutWin)
% TEP_SUBJECT_QC  Single-subject QC for the reviewer response (R1-1, R1-2).
%   - figure: ROI TEP per intensity (mean +/- trial SEM) + butterfly, removed window shaded
%   - pre-stimulus metrics per intensity: shift (-20 to -6 ms) and slope (-50 to -6 ms)
%   - residual-artefact metric: max |ROI| at 11-25 ms
% Use on data BEFORE the final 30 Hz low-pass (R1 asked for this), i.e. the _F_ERM data.

roi   = {'FC3','FC1','C3','C1'};
roiIx = find(ismember(upper({EEG.chanlocs.labels}), upper(roi)));
t     = EEG.times;
inten = nan(1, EEG.trials);                        % intensity of the time-0 TMS event
for k = 1:numel(EEG.event)
    if isfield(EEG.event, 'intensity') && ~isempty(EEG.event(k).intensity)
        inten(EEG.event(k).epoch) = EEG.event(k).intensity;
    end
end
levels = unique(inten(~isnan(inten)));

roiTr = squeeze(mean(EEG.data(roiIx, :, :), 1));   % time x trial
qc = struct('dataset', {}, 'intensity', {}, 'nTrials', {}, 'preShift_uV', {}, ...
            'preSlope_uVper10ms', {}, 'residual11to25_uV', {}, 'p30_uV', {});

fig = figure('Visible', 'off', 'Position', [100 100 1100 420]);
cols = lines(numel(levels));
subplot(1, 2, 1); hold on
for i = 1:numel(levels)
    tr = roiTr(:, inten == levels(i));
    m  = mean(tr, 2); se = std(tr, 0, 2) / sqrt(size(tr, 2));
    fill([t fliplr(t)], [m' + se' fliplr(m' - se')], cols(i, :), 'FaceAlpha', .2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t, m, 'Color', cols(i, :), 'LineWidth', 1.2, 'DisplayName', sprintf('%d%% (n=%d)', levels(i), size(tr, 2)));

    w = @(a, b) t >= a & t <= b;
    pre  = w(-20, -6); slw = w(-50, -6);
    pf   = polyfit(t(slw), m(slw)', 1);
    q.dataset = setPrefix; q.intensity = levels(i); q.nTrials = size(tr, 2);
    q.preShift_uV        = mean(m(pre)) - mean(m(w(-100, -50)));
    q.preSlope_uVper10ms = pf(1) * 10;
    q.residual11to25_uV  = max(abs(m(w(11, 25))));
    q.p30_uV             = mean(m(w(15, 50)));
    qc(end + 1) = q; %#ok<AGROW>
end
yl = ylim; patch([cutWin fliplr(cutWin)], [yl(1) yl(1) yl(2) yl(2)], [.8 .8 .8], 'EdgeColor', 'none', 'HandleVisibility', 'off');
xline(0, 'k:', 'HandleVisibility', 'off'); xlim([-100 400]); xlabel('Time from TMS pulse (ms)'); ylabel('\muV');
title(['ROI ' strjoin(roi, '/')], 'Interpreter', 'none'); legend('Location', 'best'); box off

subplot(1, 2, 2); hold on
top = levels(end);
plot(t, mean(EEG.data(:, :, inten == top), 3)', 'Color', [.6 .6 .6]);
plot(t, mean(roiTr(:, inten == top), 2), 'k', 'LineWidth', 1.5);
xline(0, 'k:'); xlim([-100 400]); xlabel('Time from TMS pulse (ms)'); ylabel('\muV');
title(sprintf('Butterfly, %d%% RMT', top)); box off
sgtitle(strrep(setPrefix, '_', '\_'));

if ~exist(fullfile(outDir, 'QC'), 'dir'), mkdir(fullfile(outDir, 'QC')); end
exportgraphics(fig, fullfile(outDir, 'QC', [setPrefix '_TEP_QC.png']), 'Resolution', 150);
close(fig);
writetable(struct2table(qc), fullfile(outDir, 'QC', [setPrefix '_TEP_QC.csv']));

for i = 1:numel(qc)
    if qc(i).residual11to25_uV > 20
        warning('%s %d%%: residual artefact %.1f uV at 11-25 ms - check ICA round 1.', setPrefix, qc(i).intensity, qc(i).residual11to25_uV);
    end
end
end
