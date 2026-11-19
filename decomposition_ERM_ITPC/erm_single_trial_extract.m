function erm_single_trial_extract()
% Single-trial ERM (ROI FC3/FC1/C3/C1, sum of IMFs 4-7) mean amplitude per epoch for the single-trial mixed model (master plan F,
% decisions 2026-10-16): P30 = 27-47 ms (ERM grand-average peak 37 ms +-10; Beck 2024 range 20-40), N100 = 120-140 ms.
% 17 participants x EO/EC, layer-1 files HHSA\timf_spm_<id>.mat (trial order = EEG epoch order of _F_ERM).
% Output: data\ERM_single_trial.csv (id, epoch, intensity label, P30, N100)
P = hhsa_setup(); spm('defaults', 'eeg'); ROOT = fileparts(P.out);
subs = [1 2 3 7 8 9 12 15 16 17 19 21 22 24 26 27 28]; rows = {};
for s = subs, for e = {'EO', 'EC'}
    id = sprintf('spTMS_%03d_Real_%s', s, e{1}); D = spm_eeg_load(fullfile(P.out, ['timf_spm_' id '.mat']));
    roi = find(ismember(D.chanlabels, {'FC3', 'FC1', 'C3', 'C1'})); t = round(D.time * 1000);   % rounded ms (float axis)
    w1 = find(t >= 27 & t <= 47); w2 = find(t >= 120 & t <= 140); cond = D.conditions;
    for k = 1:D.ntrials
        x1 = D(roi, 4:7, w1, k); x2 = D(roi, 4:7, w2, k);
        p30 = mean(squeeze(sum(mean(x1, 1), 2))); n100 = mean(squeeze(sum(mean(x2, 1), 2)));
        rows(end + 1, :) = {id, k, str2double(cond{k}), p30, n100}; %#ok<AGROW>
    end
    fprintf('%s %d trials\n', id, D.ntrials);
end, end
T = cell2table(rows, 'VariableNames', {'id', 'epoch', 'intensity', 'P30', 'N100'});
writetable(T, fullfile(ROOT, 'data', 'ERM_single_trial.csv'));
end
