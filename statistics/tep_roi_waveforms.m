function tep_roi_waveforms(outDir)
% TEP_ROI_WAVEFORMS  ROI (FC3/FC1/C3/C1) trial-average waveform per dataset x intensity from the
% conventional pipeline files (_F_TEP = 30 Hz low-pass, _F_ERM = before low-pass), for the ERM vs TEP
% overlay (Supplementary S1). No extra baseline correction (same as the ERM).
% Output: data\TEP_roi_waveforms.csv  (id, file, intensity, nTrials, t_-2000 ... t_1999)
data_dir = 'path/to/data';   % release: set this path
toolbox_dir = 'path/to/toolboxes';   % release: set this path
if nargin < 1, outDir = [data_dir '\Revision_2026\Manuscript_Figures']; end
addpath([toolbox_dir '\eeglab2021.1']); eeglab nogui;
src = [data_dir '\Revision_2026\Preproc\Final']; roi = {'FC3','FC1','C3','C1'};
F = dir(fullfile(src, '*_F_TEP.set')); rows = {}; t = [];
for f = 1:numel(F)
    id = erase(F(f).name, '_F_TEP.set');
    for v = ["TEP", "ERM"]
        EEG = pop_loadset('filename', char(id + "_F_" + v + ".set"), 'filepath', src, 'verbose', 'off');
        inten = nan(1, EEG.trials);
        for k = 1:numel(EEG.event), if ~isempty(EEG.event(k).intensity), inten(EEG.event(k).epoch) = EEG.event(k).intensity; end, end
        ri = find(ismember(upper({EEG.chanlocs.labels}), upper(roi))); t = EEG.times;
        for lv = [80 100 120]
            w = squeeze(mean(mean(EEG.data(ri, :, inten == lv), 1), 3));
            rows(end + 1, :) = [{id, char(v), lv, sum(inten == lv)}, num2cell(w)]; %#ok<AGROW>
        end
    end
    if mod(f, 10) == 0, fprintf('%d/%d\n', f, numel(F)); end
end
T = cell2table(rows, 'VariableNames', [{'id', 'file', 'intensity', 'nTrials'}, arrayfun(@(x) sprintf('t_%d', x), t, 'uni', 0)]);
writetable(T, fullfile(outDir, 'data', 'TEP_roi_waveforms.csv')); fprintf('wrote %d rows\n', height(T));
end
