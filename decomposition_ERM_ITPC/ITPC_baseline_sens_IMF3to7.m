%% ITPC_baseline_sens_IMF3to7.m  (user 2026-11-09, Part B; new file, no existing file is changed)
% ITPC baseline sensitivity analysis, IMFs 3-7 only.
% Existing analysis (structure followed here):
%   itpc_roi_export.m  : ROI ITPC (FC3/FC1/C3/C1), IMF 1-11 x -500..1000 ms x 3 intensities -> data\ITPC_roi_real.mat (X, baseline -500..-10 ms)
%   itpc_v2.py         : rule v2 = re-baseline -500..-100 ms, IMF3-7, 2 (eye) x 3 (intensity) rmANOVA with Greenhouse-Geisser correction per IMF x time,
%                        BH-FDR q = .05, one family per effect over IMF3-7 x 10..500 ms (statistics functions rm2 / bh from itpc_round15.py).
% This script re-implements those steps inline in MATLAB (no helper functions) and runs two configurations:
%   config 1 (checkpoint B2)  : baseline -500..-100 ms, tests 10..500 ms, FDR over IMF3-7 x 10..500 ms  -> expected intensity effect IMF6 27-456, IMF7 33-471 ms
%   config 2 (sensitivity B3) : baseline -500..-10 ms,  tests -100..500 ms, FDR over IMF3-7 x -100..500 ms
% Input (read only): data\ITPC_roi_real.mat, ..\Preproc\QC\QC_summary_wide.csv, data\ITPC_v2_real_anova_maps.csv (checkpoint)
% Output: data\ITPC_baseline_sens_IMF3to7.csv
clear; clc;
data_dir = 'path/to/data';   % release: set this path

%% 1. Paths and settings
MF      = [data_dir '\Revision_2026\Manuscript_Figures'];
inMat   = fullfile(MF, 'data', 'ITPC_roi_real.mat');
qcFile  = [data_dir '\Revision_2026\Preproc\QC\QC_summary_wide.csv'];
refMaps = fullfile(MF, 'data', 'ITPC_v2_real_anova_maps.csv');
outCsv  = fullfile(MF, 'data', 'ITPC_baseline_sens_IMF3to7.csv');
IMFS    = 3:7;
qFDR    = 0.05;
cfgName = {'checkpoint_v2', 'sensitivity'};
cfgBL   = [-500 -100; -500 -10];          % baseline window (ms)
cfgWin  = [-100 500; -100 500];           % analysis window (ms), as in itpc_v2.py
cfgTest = [10 500; -100 500];             % tested points = FDR family (ms)

%% 2. Participants (same rule as itpc_round15.py: processed, Real, stageD ~= exclude, both EO and EC)
Q    = readtable(qcFile, 'TextType', 'string', 'Encoding', 'UTF-8', 'VariableNamingRule', 'preserve');
keep = Q.status == "processed" & Q.stim == "Real" & Q.stageD ~= "exclude";
subs = intersect(Q.sub(keep & Q.eyes == "EO"), Q.sub(keep & Q.eyes == "EC"))';
assert(isequal(subs, [1 2 3 7 8 9 12 15 16 17 19 21 22 24 26 27 28]), 'participant list differs from the main analysis');
nSub = numel(subs);
fprintf('n = %d participants\n', nSub);

%% 3. ROI ITPC (dataset x IMF x time x intensity), rounded time axis
S   = load(inMat, 'X', 'ids', 't');
Xr  = double(S.X);
ids = cellstr(S.ids);
t   = round(double(S.t(:)'));
assert(numel(t) == size(Xr, 3));
Yall = nan(nSub, 2, 3, numel(IMFS), numel(t));            % participant x eye (EO, EC) x intensity (80, 100, 120) x IMF x time
eyes = {'EO', 'EC'};
for s = 1:nSub
    for e = 1:2
        k = find(strcmp(ids, sprintf('spTMS_%03d_Real_%s', subs(s), eyes{e})));
        assert(numel(k) == 1);
        Yall(s, e, :, :, :) = permute(squeeze(Xr(k, IMFS, :, :)), [3 1 2]);   % intensity x IMF x time
    end
end

%% 4. rmANOVA (GG) per IMF x time and BH-FDR, for each configuration
RES = struct();
for c = 1:2
    bl   = t >= cfgBL(c, 1) & t <= cfgBL(c, 2);
    tw   = t >= cfgWin(c, 1) & t <= cfgWin(c, 2);
    tt   = t(tw);
    test = tt >= cfgTest(c, 1) & tt <= cfgTest(c, 2);
    Y    = Yall(:, :, :, :, tw) - mean(Yall(:, :, :, :, bl), 5);   % baseline: subtract the mean over the baseline window
    nT   = numel(tt); nI = numel(IMFS); P = nI * nT;
    Y    = reshape(Y, nSub, 2, 3, P);                               % columns = IMF x time (IMF fastest)
    n = nSub; a = 2; b = 3;
    GM  = mean(Y, [1 2 3]); mS = mean(Y, [2 3]); mA = mean(Y, [1 3]); mB = mean(Y, [1 2]);
    mAS = mean(Y, 3); mBS = mean(Y, 2); mAB = mean(Y, 1);
    SSA   = b * n * sum((mA - GM).^2, 2);
    SSAS  = b * sum((mAS - mA - mS + GM).^2, [1 2]);
    SSB   = a * n * sum((mB - GM).^2, 3);
    SSBS  = a * sum((mBS - mB - mS + GM).^2, [1 3]);
    SSAB  = n * sum((mAB - mA - mB + GM).^2, [2 3]);
    SSABS = sum((Y - mAB - mAS - mBS + mA + mB + mS - GM).^2, [1 2 3]);
    effNames = {'eye', 'intensity', 'interaction'};
    SSeff = {SSA, SSB, SSAB}; SSerr = {SSAS, SSBS, SSABS}; dfEff = [1 2 2];
    for f = 1:3
        df  = dfEff(f); dfe = df * (n - 1);
        F   = reshape((SSeff{f} / df) ./ (SSerr{f} / dfe), 1, P);
        if df == 1
            ep = ones(1, P);
        else
            if f == 2, Z = reshape(mBS, n, 3, P); else, Z = reshape(Y(:, 1, :, :) - Y(:, 2, :, :), n, 3, P); end
            Zc = Z - mean(Z, 1);
            C  = zeros(3, 3, P);
            for i = 1:3
                for j = 1:3
                    C(i, j, :) = reshape(sum(Zc(:, i, :) .* Zc(:, j, :), 1) / (n - 1), 1, 1, P);
                end
            end
            Cd = C - mean(C, 1) - mean(C, 2) + mean(C, [1 2]);
            ep = reshape((Cd(1, 1, :) + Cd(2, 2, :) + Cd(3, 3, :)).^2 ./ ((3 - 1) * sum(Cd.^2, [1 2])), 1, P);
            ep = min(max(ep, 1 / (3 - 1)), 1);
        end
        d1 = df * ep; d2 = dfe * ep;
        p  = betainc(d2 ./ (d2 + d1 .* F), d2 / 2, d1 / 2);      % upper tail of F(d1, d2)
        F  = reshape(F, nI, nT); p = reshape(p, nI, nT); ep = reshape(ep, nI, nT);
        % BH-FDR over the tested IMF x time points
        q  = nan(nI, nT);
        pv = p(:, test); pv = pv(:); m = numel(pv);
        [ps, o] = sort(pv);
        qs = flipud(cummin(flipud(ps .* m ./ (1:m)')));
        qv = zeros(m, 1); qv(o) = min(qs, 1);
        q(:, test) = reshape(qv, nI, nnz(test));
        RES(c).(effNames{f}) = struct('F', F, 'p', p, 'eps', ep, 'q', q, 'sig', ~isnan(q) & q < qFDR);
    end
    RES(c).tt = tt; RES(c).test = test;
    fprintf('config %s: baseline %d..%d ms, window %d..%d ms, tested %d..%d ms (%d points x %d IMFs)\n', cfgName{c}, cfgBL(c, :), cfgWin(c, :), cfgTest(c, :), nnz(test), nI);
end

%% 5. Checkpoint B2: config 1 against the published v2 analysis (IMF6 27-456 ms, IMF7 33-471 ms) and data\ITPC_v2_real_anova_maps.csv
tt  = RES(1).tt;
exp67 = [27 456; 33 471];
ok  = true;
for r = 1:2
    imf = 5 + r; sig = RES(1).intensity.sig(IMFS == imf, :);
    dd  = diff([0 sig 0]); st = find(dd == 1); en = find(dd == -1) - 1;
    [~, L] = max(en - st);
    fprintf('checkpoint IMF%d intensity: longest %d-%d ms (expected %d-%d)\n', imf, tt(st(L)), tt(en(L)), exp67(r, :));
    ok = ok && tt(st(L)) == exp67(r, 1) && tt(en(L)) == exp67(r, 2);
end
R   = readtable(refMaps, 'TextType', 'string');
effNames = {'eye', 'intensity', 'interaction'};
maxdF = 0; maxdp = 0; nSigDiff = 0;
for f = 1:3
    for i = 1:numel(IMFS)
        rr = R(R.effect == effNames{f} & R.IMF == IMFS(i), :);
        assert(isequal(round(rr.time_ms'), tt));
        maxdF = max(maxdF, max(abs(rr.F' - RES(1).(effNames{f}).F(i, :))));
        maxdp = max(maxdp, max(abs(rr.pGG' - RES(1).(effNames{f}).p(i, :))));
        sigRef = ~isnan(rr.qBH') & rr.qBH' < qFDR;
        nSigDiff = nSigDiff + nnz(sigRef ~= RES(1).(effNames{f}).sig(i, :));
    end
end
fprintf('checkpoint vs ITPC_v2_real_anova_maps.csv: max |dF| = %.3g, max |dp| = %.3g, significance mismatches = %d\n', maxdF, maxdp, nSigDiff);
ok = ok && maxdF < 1e-8 && maxdp < 1e-8 && nSigDiff == 0;
if ~ok, error('STOP: checkpoint B2 failed; config 2 not reported'); end
disp('checkpoint B2 passed');

%% 6. Results per configuration, effect and IMF -> CSV
rows = {};
for c = 1:2
    tt = RES(c).tt; test = RES(c).test; pre = tt < 0;
    for f = 1:3
        for i = 1:numel(IMFS)
            sig = RES(c).(effNames{f}).sig(i, :);
            dd  = diff([0 sig 0]); st = find(dd == 1); en = find(dd == -1) - 1;
            if isempty(st)
                Ls = NaN; Le = NaN; Ln = 0; segTxt = "none";
            else
                [~, L] = max(en - st); Ls = tt(st(L)); Le = tt(en(L)); Ln = en(L) - st(L) + 1;
                segTxt = strjoin(compose("%d to %d",[tt(st)' tt(en)']), "; ");
            end
            rows(end + 1, :) = {string(cfgName{c}), sprintf('%d..%d', cfgBL(c, :)), sprintf('%d..%d', cfgTest(c, :)), ...
                sprintf('IMF3-7 x %d..%d ms (%d tests)', cfgTest(c, :), numel(IMFS) * nnz(test)), string(effNames{f}), IMFS(i), ...
                nnz(test), nnz(test & pre), nnz(sig), nnz(sig & pre), Ls, Le, Ln, segTxt}; %#ok<SAGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'config', 'baseline_ms', 'tested_ms', 'fdr_family', 'effect', 'IMF', 'n_tested_points', 'n_prestim_tested_points', ...
    'n_sig_points', 'n_sig_prestim_points', 'longest_start_ms', 'longest_end_ms', 'longest_n_points', 'all_segments_ms'});
writetable(T, outCsv);
fprintf('written: %s\n', outCsv);

%% 7. Summary of config 2 (B3)
T2 = T(T.config == "sensitivity", :);
for imf = [6 7]
    r = T2(T2.effect == "intensity" & T2.IMF == imf, :);
    fprintf('B3 IMF%d intensity: longest significant segment %d-%d ms (%d points)\n', imf, r.longest_start_ms, r.longest_end_ms, r.longest_n_points);
end
for imf = [5 6 7]
    r = T2(T2.effect == "intensity" & T2.IMF == imf, :);
    fprintf('B3 IMF%d intensity: %d significant points before 0 ms of %d pre-stimulus points tested\n', imf, r.n_sig_prestim_points, r.n_prestim_tested_points);
end
fprintf('B3 eye-state significant points (IMF3-7): %d\n', sum(T2.n_sig_points(T2.effect == "eye")));
fprintf('B3 interaction significant points (IMF3-7): %d\n', sum(T2.n_sig_points(T2.effect == "interaction")));
