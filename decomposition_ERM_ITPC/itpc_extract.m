function R = itpc_extract(id)
% ITPC_EXTRACT  Inter-trial phase coherence per intensity from the first-layer phase file
% (timph_spm_<id>.mat): ITPC = |mean over trials of exp(i*phase)| for every channel x IMF x sample,
% as spm_eeg_average with 'circularise' did in the original; then baseline 'Diff' (subtract the mean
% ITPC over the baseline window). Original baseline [-500 0] ms (file header; manuscript said -300..0);
% user 2026-10-12: [-500 -10] ms.
% Saves HHSA\itpc_<id>.mat with itpc (chan x imf x time x 3 intensities, baseline-corrected, single),
% itpc_raw for the ROI only, and the trial counts.
P = hhsa_setup(); spm('defaults', 'eeg');
D = spm_eeg_load(fullfile(P.out, ['timph_spm_' id '.mat']));
lab = D.chanlabels; t = D.time * 1000; levels = [80 100 120]; nI = size(D, 2);
bl = t >= -500 & t <= -10;
R = struct('id', id, 'time', t, 'levels', levels, 'chanlabels', {lab}, 'baseline_ms', [-500 -10], 'nTrials', zeros(1, 3));
R.itpc = nan(numel(lab), nI, numel(t), 3, 'single'); R.itpc_bl = nan(numel(lab), nI, 3, 'single');
for k = 1:3
    tr = D.indtrial(sprintf('%d', levels(k))); R.nTrials(k) = numel(tr); if isempty(tr), continue, end
    acc = zeros(numel(lab), nI, numel(t));
    for j = 1:numel(tr)                                         % accumulate unit phase vectors trial by trial (memory)
        ph = double(D(:, :, :, tr(j)));
        if j == 1 && k == 1
            fprintf('%s: phase value range [%.2f %.2f] (used as radians, as spm_eeg_average circularise does)\n', id, min(ph(:)), max(ph(:)));
        end
        acc = acc + exp(1i * ph);
    end
    itpc = abs(acc) / numel(tr);
    R.itpc_bl(:, :, k) = single(mean(itpc(:, :, bl), 3));
    R.itpc(:, :, :, k) = single(itpc - mean(itpc(:, :, bl), 3));
end
save(fullfile(P.out, ['itpc_' id '.mat']), '-struct', 'R', '-v7.3');
end
