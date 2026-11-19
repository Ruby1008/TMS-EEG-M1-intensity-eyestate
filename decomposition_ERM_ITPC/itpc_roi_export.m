function itpc_roi_export(H, outFile)
% Export ROI-averaged (FC3/FC1/C3/C1) baseline-corrected ITPC (<H>\itpc_<id>.mat, baseline -500..-10 ms) for all
% Real datasets: IMF 1-11 x -500..1000 ms x 3 intensities -> data\ITPC_roi_real.mat (v7, readable by scipy)
% 2026-10-16: optional H / outFile (AEP-ICA robustness: HHSA_AEPr -> data\ITPC_roi_real_AEPr.mat); defaults unchanged.
data_dir = 'path/to/data';   % release: set this path
if nargin < 1, H = [data_dir '\Revision_2026\Manuscript_Figures\HHSA']; end
if nargin < 2, outFile = [data_dir '\Revision_2026\Manuscript_Figures\data\ITPC_roi_real.mat']; end
F = dir(fullfile(H, 'itpc_spTMS_*_Real_*.mat')); ids = erase({F.name}, {'itpc_', '.mat'});
roi = {'FC3', 'FC1', 'C3', 'C1'}; X = []; nT = [];
for k = 1:numel(ids)
    R = load(fullfile(H, F(k).name), 'itpc', 'time', 'chanlabels', 'nTrials', 'levels');
    ri = ismember(R.chanlabels, roi); assert(sum(ri) == 4); w = R.time >= -500 & R.time <= 1000;
    x = squeeze(mean(double(R.itpc(ri, :, w, :)), 1));               % imf x time x 3
    X(k, :, :, :) = x; nT(k, :) = R.nTrials; t = R.time(w); %#ok<AGROW>
end
save(outFile, 'X', 'ids', 't', 'nT', 'roi', '-v7');
fprintf('%d datasets, size %s\n', numel(ids), mat2str(size(X)));
end
