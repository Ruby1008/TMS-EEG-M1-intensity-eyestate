function icaOverviewFigure(EEG, T, cutWin, ttl, file, decayInfo)
% ICAOVERVIEWFIGURE  One cell per component: topography, trial-average time course (-200..500 ms,
% 0 ms and the removed window marked), log power spectrum (1-100 Hz) and explained variance.
% Removed components are drawn in red with the reason; kept components in black.
% decayInfo (round 1): decay ratio/peak are printed for every component with ratio >= 15.

act = eeg_getdatact(EEG, 'component', 1:size(EEG.icaweights, 1));
n = size(act, 1); t = EEG.times; tw = t >= -200 & t <= 500;
nCol = 6; nRow = ceil(n / nCol);
cellW = 330; cellH = 200; topH = 60;
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', 'Position', [0 0 nCol * cellW, nRow * cellH + topH]);
W = nCol * cellW; H = nRow * cellH + topH;
annotation(fig, 'textbox', [0 1 - 40 / H 1 40 / H], 'String', ttl, 'EdgeColor', 'none', 'FontSize', 14, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'none');
chans = EEG.chanlocs(EEG.icachansind);
fqAll = [];
for i = 1:n
    r = ceil(i / nCol); c = i - (r - 1) * nCol;
    x0 = (c - 1) * cellW; y0 = H - topH - r * cellH;
    removed = T.removed(i);
    col = [0 0 0]; if removed, col = [0.8 0 0]; end
    lab = sprintf('IC%d  %.1f%%', i, T.varPct(i));
    if strlength(T.reason(i)) > 0, lab = sprintf('%s  [%s]', lab, T.reason(i)); end
    extra = '';
    if ~isempty(decayInfo) && ~isempty(decayInfo.ratio)
        [rmax, g] = max(decayInfo.ratio(i, :));
        if rmax >= 15, extra = sprintf('decay ratio %.0f, peak %.0f (%s)', rmax, decayInfo.peak(i, g), decayInfo.groups{g}); end
    end
    annotation(fig, 'textbox', [x0 / W, (y0 + cellH - 22) / H, cellW / W, 20 / H], 'String', lab, 'EdgeColor', 'none', ...
        'Color', col, 'FontSize', 8.5, 'FontWeight', fastif(removed, 'bold', 'normal'), 'Interpreter', 'none', 'Margin', 1);
    if ~isempty(extra)
        annotation(fig, 'textbox', [x0 / W, (y0 + cellH - 38) / H, cellW / W, 16 / H], 'String', extra, 'EdgeColor', 'none', ...
            'Color', [0.45 0.25 0], 'FontSize', 7, 'Interpreter', 'none', 'Margin', 1);
    end
    % topography
    axT = axes(fig, 'Units', 'normalized', 'Position', [(x0 + 4) / W, (y0 + 8) / H, 118 / W, 140 / H]); %#ok<LAXES>
    try
        topoplot(EEG.icawinv(:, i), chans, 'electrodes', 'off', 'gridscale', 40, 'whitebk', 'on');
    catch
        axis(axT, 'off');
    end
    % time course
    axW = axes(fig, 'Units', 'normalized', 'Position', [(x0 + 145) / W, (y0 + 92) / H, 175 / W, 62 / H]); %#ok<LAXES>
    m = mean(act(i, :, :), 3);
    hold(axW, 'on');
    yl = [min(m(tw)) max(m(tw))]; if ~(diff(yl) > 0), yl = [-1 1]; end
    yl = yl + [-0.08 0.08] * diff(yl);
    patch(axW, [cutWin fliplr(cutWin)], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], 'EdgeColor', 'none');
    plot(axW, t(tw), m(tw), 'Color', col, 'LineWidth', 0.8);
    xline(axW, 0, ':'); xlim(axW, [-200 500]); ylim(axW, yl);
    set(axW, 'FontSize', 6, 'Box', 'off', 'XTick', [-200 0 100 200 300 400 500]);
    % spectrum
    axS = axes(fig, 'Units', 'normalized', 'Position', [(x0 + 145) / W, (y0 + 14) / H, 175 / W, 58 / H]); %#ok<LAXES>
    [pxx, fq] = pwelch(double(reshape(act(i, :, :), 1, [])), EEG.srate, EEG.srate / 2, EEG.srate, EEG.srate);
    k = fq >= 1 & fq <= 100;
    plot(axS, fq(k), 10 * log10(pxx(k)), 'Color', col, 'LineWidth', 0.8);
    xlim(axS, [1 100]); set(axS, 'FontSize', 6, 'Box', 'off', 'XTick', [1 10 20 40 60 80 100]);
    if i == 1, xlabel(axW, 'ms', 'FontSize', 6); xlabel(axS, 'Hz (dB)', 'FontSize', 6); end
    fqAll = fq; %#ok<NASGU>
end
exportgraphics(fig, file, 'Resolution', 90);
close(fig);
end

function v = fastif(c, a, b)
if c, v = a; else, v = b; end
end
