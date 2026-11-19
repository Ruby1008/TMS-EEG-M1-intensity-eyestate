function erm_cluster_perm()
% ERM_CLUSTER_PERM  Whole-scalp cluster-based permutation tests of eye state on the ERM (user 2026-10-18, section 3.3; pre-specified
% analysis added for the reviewers). FieldTrip as distributed with SPM12 (<toolbox_dir>\holo_move12_tb_20240530F\holo_move12_tb\spm12\external\fieldtrip).
% Data: HHSA\erm_spTMS_<sub>_Real_<eye>.mat (erm_allch: 61 channels x 3 intensities x time; IMFs 4-7 of the first-layer CEEMDAN of _F_ERM),
% 17 main participants, 0-400 ms (time rounded to ms). Neighbours: triangulation of the 2-D projected electrode positions (waveguard64_woRef.ced).
%  (1) main effect of eye state: EO vs EC, intensities averaged; dependent-samples t, cluster alpha .05, two-tailed (correcttail 'prob'), 1000 permutations;
%  (2) eye x intensity interaction: EO - EC difference wave at 80/100/120%, dependent-samples univariate F over the three levels, cluster alpha .05
%      (F: upper tail), 1000 permutations.
% Single computational thread, no parallel pool (HHCFPC is running). rng(20261018) before each test.
% Outputs: data/ERM_cluster_perm.mat (stat structures, mask, time, labels), data/ERM_cluster_perm_clusters.csv
data_dir = 'path/to/data';   % release: set this path
toolbox_dir = 'path/to/toolboxes';   % release: set this path
P = hhsa_setup(); D = fullfile(fileparts(P.out), 'data');
ft = [toolbox_dir '\holo_move12_tb_20240530F\holo_move12_tb\spm12\external\fieldtrip'];
rmpath(genpath(P.tb)); addpath(fullfile(P.tb, 'spm12')); addpath(ft); ft_defaults; maxNumCompThreads(1);   % spm12 root: spm_bwlabel for clustering
subs = [1 2 3 7 8 9 12 15 16 17 19 21 22 24 26 27 28]; eyes = {'EO', 'EC'};
E = load(fullfile(P.out, sprintf('erm_spTMS_%03d_Real_EO.mat', subs(1)))); t = round(E.time(:)'); w = t >= 0 & t <= 400;
labs = cellstr(E.chanlabels); X = nan(numel(subs), 2, 3, numel(labs), sum(w));
for i = 1:numel(subs)
    for e = 1:2
        R = load(fullfile(P.out, sprintf('erm_spTMS_%03d_Real_%s.mat', subs(i), eyes{e})));
        assert(isequal(cellstr(R.chanlabels), labs)); assert(isequal(round(R.time(:)'), t));
        for k = 1:3, X(i, e, k, :, :) = squeeze(R.erm_allch(:, k, w)); end
    end
end
tt = t(w) / 1000;
% electrodes and neighbours
C = readtable([data_dir '\waveguard64_woRef.ced'], 'FileType', 'text', 'Delimiter', '\t');
C.labels = strtrim(C.labels); [ok, ix] = ismember(labs, C.labels); assert(all(ok));
elec = struct('label', {labs(:)}, 'elecpos', [C.X(ix) C.Y(ix) C.Z(ix)], 'unit', 'cm'); elec.chanpos = elec.elecpos;
cfgn = struct('method', 'triangulation', 'elec', elec, 'feedback', 'no'); nb = ft_prepare_neighbours(cfgn);
mk = @(Y) struct('label', {labs(:)}, 'time', tt, 'avg', Y, 'dimord', 'chan_time');
base = struct('method', 'montecarlo', 'correctm', 'cluster', 'clusteralpha', .05, 'clusterstatistic', 'maxsum', 'minnbchan', 0, ...
    'neighbours', nb, 'numrandomization', 1000, 'alpha', .05, 'latency', [0 .4], 'channel', {labs(:)});
n = numel(subs);
% (1) EO vs EC
EO = cell(1, n); EC = cell(1, n);
for i = 1:n, EO{i} = mk(squeeze(mean(X(i, 1, :, :, :), 3))); EC{i} = mk(squeeze(mean(X(i, 2, :, :, :), 3))); end
cfg = base; cfg.statistic = 'ft_statfun_depsamplesT'; cfg.tail = 0; cfg.clustertail = 0; cfg.correcttail = 'prob';
cfg.design = [1:n 1:n; ones(1, n) 2 * ones(1, n)]; cfg.uvar = 1; cfg.ivar = 2;
rng(20261018); S1 = ft_timelockstatistics(cfg, EC{:}, EO{:});        % t > 0: EC > EO
% (2) interaction: EO - EC difference wave across intensities
DF = cell(3, n);
for k = 1:3, for i = 1:n, DF{k, i} = mk(squeeze(X(i, 1, k, :, :) - X(i, 2, k, :, :))); end, end
cfg = base; cfg.statistic = 'ft_statfun_depsamplesFunivariate'; cfg.tail = 1; cfg.clustertail = 1;
cfg.design = [repmat(1:n, 1, 3); kron(1:3, ones(1, n))]; cfg.uvar = 1; cfg.ivar = 2;
rng(20261018); S2 = ft_timelockstatistics(cfg, DF{1, :}, DF{2, :}, DF{3, :});
% cluster table
rows = {};
for s = {S1, 'eye state (EC vs EO, intensities pooled)'; S2, 'eye x intensity (F on EO - EC across intensities)'}'
    S = s{1};
    for sg = {'pos', 'neg'}
        if ~isfield(S, [sg{1} 'clusters']) || isempty(S.([sg{1} 'clusters'])), continue, end
        cl = S.([sg{1} 'clusters']); lm = S.([sg{1} 'clusterslabelmat']);
        for c = 1:numel(cl)
            m = lm == c; [ci, ti] = find(m);
            rows(end + 1, :) = {s{2}, sg{1}, c, cl(c).clusterstat, cl(c).prob, 1000 * S.time(min(ti)), 1000 * S.time(max(ti)), numel(unique(ci)), strjoin(S.label(unique(ci))', ' ')}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'test', 'sign', 'cluster', 'clusterstat', 'p', 'start_ms', 'end_ms', 'n_channels', 'channels'});
writetable(T, fullfile(D, 'ERM_cluster_perm_clusters.csv'));
stat_eye = S1; stat_int = S2; labels = labs; time_ms = 1000 * tt; Xmean = squeeze(mean(X, 1)); %#ok<NASGU>
save(fullfile(D, 'ERM_cluster_perm.mat'), 'stat_eye', 'stat_int', 'labels', 'time_ms', 'Xmean', '-v7');
disp(T(T.p < .2, :));
end
