function R = erm_extract(id, varargin)
% ERM_EXTRACT  Event-related mode (ERM) TEP from the first-layer CEEMDAN file of one dataset,
% exactly as spTMS_ERM.m did: ROI FC3/FC1/C3/C1 (by label), IMFs 4-7, mean over ROI then over
% trials of each intensity, then sum the IMFs; no baseline correction; 20-ms sliding mean
% (21 samples, 10-ms step, 0-400 ms) for the sliding-window ANOVA.
% Also returns the all-channel ERM per intensity (for the topography inset; the original script
% used EO/100% data for all EC/120% arrays - fixed here, user 2026-10-12) and the pre-stimulus
% metrics with the Phase A definition (ref -300..-100, shift -20..-5, slope -60..-5, post 15..50 ms).
%   R = erm_extract('spTMS_001_Real_EO')   -> also saves HHSA\erm_<id>.mat
p = inputParser; p.addParameter('imfs', 4:7); p.addParameter('roi', {'FC3','FC1','C3','C1'}); p.parse(varargin{:});
P = hhsa_setup(); spm('defaults', 'eeg');
D = spm_eeg_load(fullfile(P.out, ['timf_spm_' id '.mat']));
lab = D.chanlabels; roi = find(ismember(lab, p.Results.roi)); assert(numel(roi) == 4, '%s: ROI channels missing', id);
t = D.time * 1000; imfs = p.Results.imfs;
levels = [80 100 120]; nT = zeros(1, 3);
R = struct('id', id, 'time', t, 'levels', levels, 'imfs', imfs, 'roi', {lab(roi)}, 'chanlabels', {lab});
R.erm_roi = nan(3, numel(t)); R.erm_allch = nan(numel(lab), 3, numel(t)); R.erm_roi_imf = nan(3, numel(imfs), numel(t));
for k = 1:3
    tr = D.indtrial(sprintf('%d', levels(k))); nT(k) = numel(tr); if isempty(tr), continue, end
    x = D(roi, imfs, :, tr);                                   % roi x imf x time x trial
    m = mean(mean(x, 1), 4);                                    % 1 x imf x time
    R.erm_roi_imf(k, :, :) = m;
    R.erm_roi(k, :) = squeeze(sum(m, 2))';
    a = mean(D(:, imfs, :, tr), 4);                             % chan x imf x time  (trial average)
    R.erm_allch(:, k, :) = squeeze(sum(a, 2));
end
R.nTrials = nT;
% 20-ms sliding mean, 10-ms step, 0..400 ms (spTMS_ERM.m: mean(x(t:t+20)), t = 0:10:390)
starts = 0:10:390; R.win_start = starts; R.erm_roi_win = nan(3, numel(starts));
for w = 1:numel(starts)
    i0 = find(t >= starts(w), 1); R.erm_roi_win(:, w) = mean(R.erm_roi(:, i0:i0 + 20), 2);
end
% pre-stimulus metrics on the ERM (Phase A definition, as in Stats_R1-1)
R.prestim = struct('intensity', num2cell(levels), 'shift', [], 'slope', [], 'post', []);
for k = 1:3
    w = R.erm_roi(k, :); w = w - mean(w(t >= -300 & t <= -100));
    s = t >= -60 & t <= -5; pf = polyfit(t(s), w(s), 1);
    R.prestim(k).shift = mean(w(t >= -20 & t <= -5)); R.prestim(k).slope = pf(1) * 10; R.prestim(k).post = mean(w(t >= 15 & t <= 50));
end
save(fullfile(P.out, ['erm_' id '.mat']), '-struct', 'R');
end
