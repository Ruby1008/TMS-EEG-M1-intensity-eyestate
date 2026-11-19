% CHECK_TRIALS_VS_2023  Are the trials of the new (.cnt, pulse-locked) preprocessing the same as the 2023 _e epochs?
% READ-ONLY on the 2023 files. Each TMS pulse has a unique artefact size (max |diff| at FC3 around the pulse);
% both versions come from the same raw recording, so the same trial gives the same value. The two pulse-size
% sequences are matched one to one (tolerance 0.5 uV), and the intensity labels of matched trials are compared.
% Output: Preproc\trials_vs_2023.csv (one row per dataset) and trials_vs_2023_detail.csv (one row per trial).
data_dir = 'path/to/data';   % release: set this path
toolbox_dir = 'path/to/toolboxes';   % release: set this path
addpath([toolbox_dir '\eeglab2021.1']); eeglab nogui; close all
src = data_dir;
out = [data_dir '\Revision_2026\Preproc'];
if ~exist('listName', 'var'), listName = 'datalist_Real.csv'; end     % set listName = 'datalist_Sham.csv' before running for Sham
if ~exist('pulseThr', 'var'), pulseThr = 1000; end                     % 300 for Sham (findpulse rate used for Sham)
if ~exist('outTag', 'var'), outTag = ''; end                          % '_Sham' -> trials_vs_2023_Sham.csv / _detail_Sham.csv
L = readtable(fullfile(out, listName), 'TextType', 'string', 'Delimiter', ',');
L = L(L.exclude == 0, :);
R = []; Dt = [];
for i = 1:height(L)
    id = char(L.id(i)); r = struct('id', string(id), 'e2023_trials', NaN, 'new_pulses', NaN, 'matched', NaN, ...
        'only2023', NaN, 'onlyNew', NaN, 'sameOrder', NaN, 'labelAgree', NaN, 'labelCompared', NaN, 'note', "");
    P = load(fullfile(out, 'Logs', [id '_pairInfo.mat'])); pi = P.pairInfo;
    newSize = pi.pulseSize(:)'; newInt = pi.intensity(:)'; r.new_pulses = numel(newSize);
    ef = fullfile(src, [id '_e.set']);
    if ~exist(ef, 'file'), r.note = "no 2023 _e file"; R = [R r]; continue, end %#ok<AGROW>
    E = pop_loadset('filename', [id '_e.set'], 'filepath', src);
    t = E.times; ch = strcmpi({E.chanlocs.labels}, 'FC3'); r.e2023_trials = E.trials;
    ev = string(L.evCodes(i)); codes = num2cell(char(ev));
    oldSize = nan(1, E.trials); oldInt = nan(1, E.trials);
    cl = cell(1, E.trials);
    for k = 1:numel(E.event)
        ty = num2str(E.event(k).type); j = find(strcmp(codes, ty));
        if ~isempty(j), lk = mod(E.event(k).latency - 1, E.pnts) + 1; cl{E.event(k).epoch}(end + 1, :) = [t(round(lk)) j]; end
    end
    for tr = 1:E.trials
        x = double(E.data(ch, :, tr)); dx = abs(diff(x));
        w = find(t(1:end - 1) >= -20 & t(1:end - 1) <= 150); k = find(dx(w) > pulseThr, 1);
        if isempty(k), continue, end
        s = w(k) + 1; ww = max(1, s - 2):min(E.pnts, s + 5);
        oldSize(tr) = max(abs(diff(x(ww))));
        if ~isempty(cl{tr})
            dl = t(s) - cl{tr}(:, 1); dl(dl < 0 | dl > 100) = Inf; [m, j] = min(dl);
            if isfinite(m), oldInt(tr) = 60 + 20 * cl{tr}(j, 2); end     % 1/2/3 -> 80/100/120
        end
    end
    % one-to-one matching on pulse size
    used = false(size(newSize)); match = nan(1, E.trials);
    for tr = 1:E.trials
        if isnan(oldSize(tr)), continue, end
        d = abs(newSize - oldSize(tr)); d(used) = Inf; [m, j] = min(d);
        if m <= 0.5, match(tr) = j; used(j) = true; end
    end
    mOK = ~isnan(match);
    r.matched = sum(mOK); r.only2023 = sum(~mOK); r.onlyNew = sum(~used);
    r.sameOrder = all(diff(match(mOK)) > 0);
    cmp = mOK & ~isnan(oldInt); r.labelCompared = sum(cmp);
    r.labelAgree = sum(oldInt(cmp) == newInt(match(cmp)));
    if r.labelCompared == 0, r.note = "2023 _e has no trigger codes (labels not compared)"; end
    R = [R r]; %#ok<AGROW>
    for tr = 1:E.trials
        nI = NaN; if mOK(tr), nI = newInt(match(tr)); end
        Dt = [Dt struct('id', string(id), 'trial2023', tr, 'newTrial', match(tr), 'pulseSize2023', oldSize(tr), 'int2023', oldInt(tr), 'intNew', nI)]; %#ok<AGROW>
    end
    fprintf('%s: 2023 %d | new %d | matched %d | only2023 %d | onlyNew %d | order %d | labels %d/%d\n', id, r.e2023_trials, ...
        r.new_pulses, r.matched, r.only2023, r.onlyNew, r.sameOrder, r.labelAgree, r.labelCompared);
end
writetable(struct2table(R), fullfile(out, ['trials_vs_2023' outTag '.csv']));
writetable(struct2table(Dt), fullfile(out, ['trials_vs_2023_detail' outTag '.csv']));
