function imf_frequency_distribution(which, varargin)
% IMF_FREQUENCY_DISTRIBUTION  Fig 4: distribution of the instantaneous frequency of each IMF (first-layer
% timfm_ files = instantaneous frequency per channel x IMF x sample x trial), averaged over all electrodes,
% trials and participants. The original was produced in the toolbox GUI (FD*.mat); here: probability
% histogram on a log2 frequency axis (0.25-512 Hz, 1/8-octave bins), amplitude-unweighted, per IMF;
% dashed line = mean of the distribution, label = peak (mode) in Hz, as in the manuscript figure.
%   imf_frequency_distribution('Real')      datasets: all Real with timfm_ present (set 'ids' to restrict)
% Output: Fig4_new.png/.fig, data\Fig4_IMF_frequency_distribution.csv (bin x IMF probability), data\Fig4_IMF_peak_mean_Hz.csv
p = inputParser; p.addParameter('ids', {}); p.addParameter('imfs', 1:10); p.addParameter('tag', ''); p.parse(varargin{:});
P = hhsa_setup(); spm('defaults', 'eeg'); ROOT = fileparts(P.out);
if isempty(p.Results.ids), F = dir(fullfile(P.out, ['timfm_spm_spTMS_*_' which '_*.mat'])); ids = erase({F.name}, {'timfm_spm_', '.mat'}); else, ids = p.Results.ids; end
edges = 2 .^ (-2:1/8:9); ctr = sqrt(edges(1:end - 1) .* edges(2:end)); imfs = p.Results.imfs;
H = zeros(numel(ctr), numel(imfs)); M = zeros(numel(ids), numel(imfs)); n = 0;
for k = 1:numel(ids)
    D = spm_eeg_load(fullfile(P.out, ['timfm_spm_' ids{k} '.mat'])); ti = D.indsample(-1):D.indsample(1);   % -1..1 s (edges excluded)
    Hk = zeros(numel(ctr), numel(imfs));
    for tr = 1:D.ntrials
        x = D(:, imfs, ti, tr);
        for j = 1:numel(imfs)
            v = x(:, j, :); v = v(v > 0 & isfinite(v)); Hk(:, j) = Hk(:, j) + histcounts(v, edges)';
        end
    end
    Hk = Hk ./ sum(Hk, 1); H = H + Hk; n = n + 1; M(k, :) = sum(Hk .* log2(ctr)', 1);
    fprintf('%d/%d %s\n', k, numel(ids), ids{k});
end
H = H / n; pk = ctr(arrayfun(@(j) find(H(:, j) == max(H(:, j)), 1), 1:numel(imfs))); mn = 2 .^ mean(M, 1);
tag = p.Results.tag; if isempty(tag), tag = which; end
writetable(array2table([ctr' H], 'VariableNames', [{'freq_Hz'}, arrayfun(@(j) sprintf('IMF%d', j), imfs, 'uni', 0)]), fullfile(ROOT, 'data', ['Fig4_IMF_frequency_distribution_' tag '.csv']));
writetable(table(imfs', pk', mn', 'VariableNames', {'IMF', 'peak_Hz', 'mean_Hz'}), fullfile(ROOT, 'data', ['Fig4_IMF_peak_mean_Hz_' tag '.csv']));
% figure in the style of the manuscript Fig 4
cols = [0 114 189; 126 47 142; 119 172 48; 77 190 238; 0 0 255; 0 128 0; 204 0 204; 190 190 0; 80 80 80; 0 102 204] / 255;
fig = figure('Visible', 'off', 'Position', [0 0 900 650], 'Color', 'w'); hold on
for j = 1:numel(imfs)
    c = cols(mod(j - 1, size(cols, 1)) + 1, :); plot(log2(ctr), H(:, j), 'Color', c, 'LineWidth', 3, 'DisplayName', sprintf('IMF%d', imfs(j)));
    xline(log2(mn(j)), '--', 'Color', c, 'LineWidth', 2, 'HandleVisibility', 'off');
    text(log2(pk(j)), max(H(:, j)) + .01, sprintf('%.2f', pk(j)), 'Color', c, 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'center');
end
xt = [0.25 0.5 1 2 4 8 16 32 64 128 256]; set(gca, 'XTick', log2(xt), 'XTickLabel', arrayfun(@num2str, xt, 'uni', 0), 'FontSize', 16, 'LineWidth', 1.5, 'FontName', 'Times New Roman');
xlim(log2([0.2 512])); ylim([0 0.5]); xlabel('Frequency', 'FontSize', 20); ylabel('Probability', 'FontSize', 20); title('IMF Distribution', 'FontSize', 22, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 14); box on
fn = fullfile(ROOT, ['Fig4_new_' tag]); print(fig, [fn '.png'], '-dpng', '-r600'); print(fig, [fn '.tif'], '-dtiff', '-r600'); savefig(fig, [fn '.fig']); close(fig)
fprintf('IMF peak Hz: %s\nIMF mean Hz: %s\n', mat2str(round(pk, 2)), mat2str(round(mn, 2)));
end
