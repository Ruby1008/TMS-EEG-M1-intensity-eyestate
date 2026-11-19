function extract_prestim_phaseA_def(outDir)
% EXTRACT_PRESTIM_PHASEA_DEF  Pre-stimulus metrics with the Phase A definition (R1-1; user 2026-10-12),
% so the new pipeline can be compared with the original manuscript (phaseA_diagnostics.m, prestim_metrics):
%   ROI FC3/FC1/C3/C1 trial average per intensity, minus mean(-300..-100 ms);
%   shift = mean(-20..-5 ms), slope = -60..-5 ms (uV/10 ms), post = mean(15..50 ms)
% Run on Final\<id>_F_TEP.set (main) and _F_ERM.set (no 30 Hz low-pass).
% Output: Stats_R1-1\prestim_phaseA_def.csv (one row per dataset x file x intensity)
data_dir = 'path/to/data';   % release: set this path
toolbox_dir = 'path/to/toolboxes';   % release: set this path
if nargin < 1, outDir = [data_dir '\Revision_2026\Preproc']; end
addpath([toolbox_dir '\eeglab2021.1']); eeglab nogui
roi = {'FC3', 'FC1', 'C3', 'C1'};
F = dir(fullfile(outDir, 'Final', '*_F_TEP.set'));
rows = {};
for f = 1:numel(F)
    id = erase(F(f).name, '_F_TEP.set');
    for v = ["TEP", "ERM"]
        EEG = pop_loadset('filename', char(id + "_F_" + v + ".set"), 'filepath', fullfile(outDir, 'Final'), 'verbose', 'off');
        inten = nan(1, EEG.trials);
        for k = 1:numel(EEG.event)
            if isfield(EEG.event, 'intensity') && ~isempty(EEG.event(k).intensity), inten(EEG.event(k).epoch) = EEG.event(k).intensity; end
        end
        ri = find(ismember(upper({EEG.chanlocs.labels}), upper(roi)));
        t = EEG.times;
        for lv = unique(inten(~isnan(inten)))
            w = squeeze(mean(mean(EEG.data(ri, :, inten == lv), 1), 3));
            w = w - mean(w(t >= -300 & t <= -100));
            s = t >= -60 & t <= -5; p = polyfit(t(s), w(s), 1);
            rows(end + 1, :) = {id, char(v), lv, sum(inten == lv), numel(ri), mean(w(t >= -20 & t <= -5)), p(1) * 10, mean(w(t >= 15 & t <= 50))}; %#ok<AGROW>
        end
    end
    fprintf('%d/%d %s\n', f, numel(F), id);
end
T = cell2table(rows, 'VariableNames', {'id', 'file', 'intensity', 'nTrials', 'nRoiChan', 'shift_m20_m5', 'slope_m60_m5_uVper10ms', 'post_15_50'});
if ~exist(fullfile(outDir, 'Stats_R1-1'), 'dir'), mkdir(fullfile(outDir, 'Stats_R1-1')); end
writetable(T, fullfile(outDir, 'Stats_R1-1', 'prestim_phaseA_def.csv'));
fprintf('wrote %d rows\n', height(T));
end
