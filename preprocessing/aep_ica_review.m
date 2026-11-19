function aep_ica_review(subs)
% AEP_ICA_REVIEW  AEP-ICA robustness analysis, step 1 (ICA + review figures; NO removal).
% Replicates the original procedure in <data_dir>/spTMS_AEP_20241128 (read-only reference):
%   four blocks concatenated in the order Real EC, Real EO, Sham EC, Sham EO (merge.m)  ->
%   extended Infomax (pop_runica 'extended',1) with PCA to 40 dimensions, computed on -0.5..0.6 s and applied
%   to the whole epoch (AEP_ICs_Loop*.m; confirmed for 18/23 original *_All4RS_mi40.set files by the PCA subspace
%   of the -0.5..0.6 s data, Revision_2026\AEP_inventory_2026-10-14.html)  ->  pop_iclabel(EEG,'default')  ->
%   visual selection (here: review page, user ticks; removal is a separate step).
% New-pipeline input (user 2026-10-14): Preproc\Intermediate\<id>_eririffiiiRe.set = after average reference,
% before baseline and before the 30 Hz low-pass; one ICA per participant; the same removal will later give
% _F_ERM_AEPr and _F_TEP_AEPr. Only the blocks that exist are concatenated. rng fixed (runica rndreset off).
% Output: Revision_2026\Preproc_AEPr\ICA\<sub>_AEP_ICA.mat, Preproc_AEPr\Review\figs\<sub>_IC##.png (300 dpi),
%         Preproc_AEPr\Review\<sub>_review.json
data_dir = 'path/to/data';   % release: set this path
toolbox_dir = 'path/to/toolboxes';   % release: set this path
if nargin < 1, subs = [1:4 6:9 12:13 15:19 21:24 26:28]; end
addpath([toolbox_dir '\eeglab2021.1']); eeglab nogui;
pre = [data_dir '\Revision_2026\Preproc\Intermediate'];
out = [data_dir '\Revision_2026\Preproc_AEPr'];
icaDir = fullfile(out, 'ICA'); figDir = fullfile(out, 'Review', 'figs');
if ~exist(icaDir, 'dir'), mkdir(icaDir); end, if ~exist(figDir, 'dir'), mkdir(figDir); end
blocks = {'Real_EC', 'Real_EO', 'Sham_EC', 'Sham_EO'};
roi = {'FC3', 'FC1', 'C3', 'C1'};
bcol = [57 106 106; 197 85 5; 150 150 150; 0 0 0] / 255;          % Real EC, Real EO, Sham EC, Sham EO
set(groot, 'defaultAxesFontName', 'Arial', 'defaultTextFontName', 'Arial');
for s = subs
    sid = sprintf('%03d', s); f = fullfile(icaDir, ['spTMS_' sid '_AEP_ICA.mat']);
    if exist(f, 'file'), fprintf('%s exists, skip\n', sid); continue, end
    t0 = tic; ALL = []; nb = []; have = {};
    for b = 1:4
        fn = sprintf('spTMS_%s_%s_eririffiiiRe.set', sid, blocks{b});
        if ~exist(fullfile(pre, fn), 'file'), continue, end
        E = pop_loadset('filename', fn, 'filepath', pre, 'verbose', 'off');
        if isempty(ALL), ALL = E; else, ALL = pop_mergeset(ALL, E, 0); end              % pairwise merge (block structs differ in fields)
        nb(end + 1) = E.trials; have{end + 1} = blocks{b}; %#ok<AGROW>
    end
    if numel(nb) < 2, fprintf('%s: fewer than 2 blocks, skip\n', sid); continue, end
    EEG = ALL; assert(EEG.trials == sum(nb), 'trial count mismatch after merge');
    blk = repelem(1:numel(nb), nb);                                   % block index per epoch
    Et = pop_select(EEG, 'time', [-0.5 0.6]);
    rng(20261014);
    Et = pop_runica(Et, 'icatype', 'runica', 'extended', 1, 'pca', 40, 'rndreset', 'no', 'interrupt', 'off');
    EEG.icaweights = Et.icaweights; EEG.icasphere = Et.icasphere; EEG.icawinv = Et.icawinv; EEG.icachansind = Et.icachansind; EEG.icaact = [];
    EEG = eeg_checkset(EEG, 'ica');
    EEG = pop_iclabel(EEG, 'default');
    cls = EEG.etc.ic_classification.ICLabel.classes; P = EEG.etc.ic_classification.ICLabel.classifications;
    act = eeg_getdatact(EEG, 'component', 1:size(EEG.icaweights, 1));
    t = EEG.times; lab = {EEG.chanlocs.labels}; iCz = find(strcmp(lab, 'Cz')); iR = find(ismember(lab, roi));
    pv = zeros(1, size(act, 1)); X = reshape(EEG.data, EEG.nbchan, []);
    for c = 1:size(act, 1), proj = EEG.icawinv(:, c) * reshape(act(c, :, :), 1, []); pv(c) = 100 * (1 - var(X(:) - proj(:)) / var(X(:))); end
    info = struct('sub', sid, 'blocks', {have}, 'nTrials', nb, 'nIC', size(act, 1), 'iclabelClasses', {cls}, 'ics', []);
    for c = 1:size(act, 1)
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 1500 330]);
        % topography
        axes('Position', [0.0 0.06 0.12 0.80]); topoplot(EEG.icawinv(:, c), EEG.chanlocs(EEG.icachansind), 'electrodes', 'off'); %#ok<LAXES>
        % ICLabel bar
        ax = axes('Position', [0.155 0.18 0.09 0.66]); barh(ax, P(c, :), 'FaceColor', [0.4 0.4 0.4]); %#ok<LAXES>
        set(ax, 'YTick', 1:numel(cls), 'YTickLabel', cls, 'FontSize', 6, 'XLim', [0 1], 'YDir', 'reverse', 'Box', 'off'); title(ax, 'ICLabel', 'FontSize', 8);
        % ERP image (all epochs, blocks stacked), -200..500 ms
        w = t >= -200 & t <= 500; A = squeeze(act(c, w, :))'; A = conv2(A, ones(5, 1) / 5, 'same');
        ax = axes('Position', [0.29 0.18 0.17 0.66]); imagesc(ax, t(w), 1:size(A, 1), A); hold(ax, 'on'); %#ok<LAXES>
        cl = prctile(abs(A(:)), 98); clim(ax, [-cl cl]); colormap(ax, parula);
        for b = 1:numel(nb) - 1, yline(ax, sum(nb(1:b)) + 0.5, 'w-', 'LineWidth', 1); end
        xline(ax, 0, 'w:'); set(ax, 'YTick', cumsum(nb) - nb / 2, 'YTickLabel', strrep(have, '_', ' '), 'FontSize', 6);
        title(ax, 'ERP image（epoch，5-trial smooth）', 'FontSize', 8); xlabel(ax, 'ms', 'FontSize', 7);
        % spectrum
        ax = axes('Position', [0.50 0.18 0.11 0.66]); %#ok<LAXES>
        [pxx, fq] = pwelch(double(reshape(act(c, :, :), 1, [])), EEG.srate, EEG.srate / 2, EEG.srate, EEG.srate); k = fq >= 1 & fq <= 80;
        plot(ax, fq(k), 10 * log10(pxx(k)), 'k'); xlim(ax, [1 80]); title(ax, '頻譜（dB）', 'FontSize', 8); set(ax, 'FontSize', 7, 'Box', 'off');
        % back-projection to Cz and ROI, block averages
        for q = 1:2
            ax = axes('Position', [0.65 + (q - 1) * 0.175 0.18 0.155 0.66]); hold(ax, 'on'); %#ok<LAXES>
            ch = iCz; ttl = 'Cz'; if q == 2, ch = iR; ttl = 'ROI FC3/FC1/C3/C1'; end
            for b = 1:numel(nb)
                a = mean(act(c, :, blk == b), 3); y = mean(EEG.icawinv(ch, c), 1) * a;
                plot(ax, t(w), y(w), 'Color', bcol(strcmp(blocks, have{b}), :), 'LineWidth', 1.1, 'DisplayName', strrep(have{b}, '_', ' '));
            end
            xline(ax, 0, ':', 'HandleVisibility', 'off'); for L = [50 100 200], xline(ax, L, ':', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off'); end
            xlim(ax, [-200 500]); title(ax, ['反投影到 ' ttl '（µV）'], 'FontSize', 8); set(ax, 'FontSize', 7, 'Box', 'off');
            if q == 2, legend(ax, 'FontSize', 6, 'Location', 'eastoutside', 'Box', 'off'); end
        end
        [pm, im] = max(P(c, :));
        annotation(fig, 'textbox', [0 0.88 1 0.12], 'String', sprintf('spTMS %s  IC%d   ICLabel: %s %.2f   變異 %.1f%%', sid, c, cls{im}, pm, pv(c)), ...
            'EdgeColor', 'none', 'FontSize', 10, 'Interpreter', 'none');
        png = sprintf('spTMS_%s_IC%02d.png', sid, c);
        exportgraphics(fig, fullfile(figDir, png), 'Resolution', 300); close(fig);
        info.ics = [info.ics struct('ic', c, 'iclabel', {cls(im)}, 'iclabelP', P(c, :), 'varPct', pv(c), 'png', png)];
    end
    W = struct('icaweights', EEG.icaweights, 'icasphere', EEG.icasphere, 'icawinv', EEG.icawinv, 'icachansind', EEG.icachansind, ...
        'iclabel', P, 'iclabelClasses', {cls}, 'blocks', {have}, 'nTrials', nb, 'icaWindow_s', [-0.5 0.6], 'seed', 20261014, 'chanlabels', {lab});
    save(f, '-struct', 'W');
    fid = fopen(fullfile(out, 'Review', ['spTMS_' sid '_review.json']), 'w', 'n', 'UTF-8'); fwrite(fid, jsonencode(info), 'char'); fclose(fid);
    fprintf('%s: %d blocks (%s), %d ICs, %.1f min\n', sid, numel(nb), strjoin(have, ','), size(act, 1), toc(t0) / 60);
    close all force
end
end
