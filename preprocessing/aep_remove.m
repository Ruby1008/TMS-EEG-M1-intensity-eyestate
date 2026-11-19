function aep_remove(subs)
% AEP_REMOVE  AEP-ICA robustness, step 2: remove the components ticked by the user and write _F_ERM_AEPr / _F_TEP_AEPr.
% Selection: Revision_2026\Preproc_AEPr\Review\by_subject\aep_selection.csv (sub, component, action).
% User decisions 2026-10-15 (HANDOFF section 15): only 'remove' is removed ('unsure' = keep); 001 was not rated = no
% removal; the large-variance non-Brain components ticked 'remove' are removed as ticked.
% ICA: Preproc_AEPr\ICA\spTMS_<sub>_AEP_ICA.mat (aep_ica_review.m; one ICA on the concatenated blocks). Removal is applied
% block by block with the stored weights (identical to pop_subcomp on the concatenated data, which is a per-sample
% linear operation), on Preproc\Intermediate\<id>_eririffiiiRe.set (after average reference, before baseline).
% Then exactly as preprocess_tep.m: baseline -100..-6 ms -> _F_ERM_AEPr; 30 Hz 4th-order zero-phase Butterworth
% low-pass (pop_tesa_filtbutter), no second baseline -> _F_TEP_AEPr; event types '80'/'100'/'120'.
% Outputs: Preproc_AEPr\Final\<id>_F_ERM_AEPr.set / _F_TEP_AEPr.set; Preproc_AEPr\AEP_removal_log.csv;
%          Preproc_AEPr\Review\removal_check\spTMS_<sub>.png (ROI before/after per block, 150 dpi).
data_dir = 'path/to/data';   % release: set this path
toolbox_dir = 'path/to/toolboxes';   % release: set this path
REV = [data_dir '\Revision_2026'];
pre = fullfile(REV, 'Preproc', 'Intermediate'); fin0 = fullfile(REV, 'Preproc', 'Final');
A = fullfile(REV, 'Preproc_AEPr'); outF = fullfile(A, 'Final'); chk = fullfile(A, 'Review', 'removal_check');
if ~exist(outF, 'dir'), mkdir(outF); end, if ~exist(chk, 'dir'), mkdir(chk); end
addpath([toolbox_dir '\eeglab2021.1']); eeglab nogui; close all
T = readtable(fullfile(A, 'Review', 'by_subject', 'aep_selection.csv'), 'CommentStyle', '#', 'Delimiter', ',', 'TextType', 'string');
F = dir(fullfile(A, 'ICA', 'spTMS_*_AEP_ICA.mat'));
all_subs = cellfun(@(f) str2double(f(7:9)), {F.name});
if nargin < 1, subs = all_subs; end
roi = {'FC3', 'FC1', 'C3', 'C1'}; lev = [80 100 120]; cols = [0.2 0.6 0.2; 0.9 0.5 0; 0.8 0 0];
L = {};
for s = subs
    sid = sprintf('%03d', s); W = load(fullfile(A, 'ICA', ['spTMS_' sid '_AEP_ICA.mat']));
    r = T(T.sub == s, :); rm = sort(r.component(r.action == "remove"))';
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 400 * numel(W.blocks) 330]);
    for b = 1:numel(W.blocks)
        id = sprintf('spTMS_%s_%s', sid, W.blocks{b});
        EEG = pop_loadset('filename', [id '_eririffiiiRe.set'], 'filepath', pre, 'verbose', 'off');
        assert(isequal({EEG.chanlocs.labels}, W.chanlabels), '%s: channel labels differ from the ICA', id);
        assert(EEG.trials == W.nTrials(b), '%s: trial count differs from the ICA', id);
        if ~isempty(rm)
            EEG.icaweights = W.icaweights; EEG.icasphere = W.icasphere; EEG.icawinv = W.icawinv; EEG.icachansind = W.icachansind; EEG.icaact = [];
            EEG = eeg_checkset(EEG, 'ica');
            EEG = pop_subcomp(EEG, rm, 0);
        end
        EEG.icaweights = []; EEG.icasphere = []; EEG.icawinv = []; EEG.icachansind = []; EEG.icaact = [];
        EEG = eeg_checkset(EEG);
        EEG = pop_rmbase(EEG, [-100 -6]);
        ERM = EEG; TEP = pop_tesa_filtbutter(EEG, [], 30, 4, 'lowpass');
        for v = 1:2
            if v == 1, E = ERM; suf = '_F_ERM_AEPr'; else, E = TEP; suf = '_F_TEP_AEPr'; end
            E.event = E.event(strcmp(cellfun(@num2str, {E.event.type}, 'uni', 0), 'TMS'));
            for i = 1:numel(E.event), E.event(i).type = num2str(E.event(i).intensity); end
            E = eeg_checkset(E, 'eventconsistency'); E.setname = [id suf];
            pop_saveset(E, 'filename', [id suf '.set'], 'filepath', outF, 'savemode', 'twofiles');
        end
        % check against the main-pipeline output (no AEP removal)
        O = pop_loadset('filename', [id '_F_ERM.set'], 'filepath', fin0, 'verbose', 'off');
        Ot = pop_loadset('filename', [id '_F_TEP.set'], 'filepath', fin0, 'verbose', 'off');
        iR = ismember({ERM.chanlocs.labels}, roi); t = ERM.times; w = t >= -100 & t <= 400;
        inten = nan(1, ERM.trials); for k = 1:numel(ERM.event), if ~isempty(ERM.event(k).intensity), inten(ERM.event(k).epoch) = ERM.event(k).intensity; end, end
        d = double(ERM.data) - double(O.data); dt = double(TEP.data) - double(Ot.data);
        ax = subplot(1, numel(W.blocks), b); hold(ax, 'on');
        for k = 1:3
            plot(ax, t(w), mean(mean(O.data(iR, w, inten == lev(k)), 3), 1), '-', 'Color', cols(k, :), 'LineWidth', 1);
            plot(ax, t(w), mean(mean(ERM.data(iR, w, inten == lev(k)), 3), 1), '--', 'Color', cols(k, :), 'LineWidth', 1);
        end
        xline(ax, 0, ':'); title(ax, sprintf('%s（實線前、虛線後）', strrep(W.blocks{b}, '_', ' ')), 'FontSize', 8); xlabel(ax, 'ms'); ylabel(ax, 'ROI µV');
        L(end + 1, :) = {s, W.blocks{b}, EEG.trials, numel(rm), strjoin(string(rm), ' '), max(abs(d(:))), max(abs(dt(:))), ...
                         max(abs(mean(mean(d(iR, w, :), 3), 1)))}; %#ok<AGROW>
        fprintf('%s: removed %d; max|after-before| all %.3g uV, TEP %.3g\n', id, numel(rm), max(abs(d(:))), max(abs(dt(:))));
    end
    sgtitle(fig, sprintf('spTMS %s：AEP 移除 %d 個（%s），ROI FC3/FC1/C3/C1，綠 80／橙 100／紅 120%%', sid, numel(rm), strjoin(string(rm), ',')), 'FontSize', 9);
    exportgraphics(fig, fullfile(chk, ['spTMS_' sid '.png']), 'Resolution', 150); close(fig);
end
C = cell2table(L, 'VariableNames', {'sub', 'block', 'trials', 'n_removed', 'removed_ICs', 'max_abs_diff_ERM_uV', 'max_abs_diff_TEP_uV', 'max_abs_diff_ROIavg_uV'});
f = fullfile(A, 'AEP_removal_log.csv');
if exist(f, 'file') && nargin >= 1, C = [readtable(f, 'TextType', 'string'); C]; end %#ok<AGROW>
writetable(C, f);
end
