function status = run_preproc_batch(which, varargin)
% RUN_PREPROC_BATCH  Run preprocess_tep over the data list.
%   run_preproc_batch('pilot')            the 3 pilot datasets (column 'pilot' not empty)
%   run_preproc_batch('all')              every dataset with exclude == 0
%   run_preproc_batch({'spTMS_001_Real_EC'})
% Options: 'workers' (0 = serial), 'force' (re-run even if up to date)
% A dataset is re-run only when its final files are missing, or its rows in ica_overrides.csv
% changed since the last run (README 6b). ICA weights are reused, so component numbers stay the same.
data_dir = 'path/to/data';   % release: set this path

p = inputParser;
p.addParameter('workers', 0);
p.addParameter('force', false);
p.addParameter('outDir', [data_dir '\Revision_2026\Preproc']);
p.addParameter('extra', {});              % extra name/value pairs for preprocess_tep (test variants only)
p.addParameter('fromE', false);
p.addParameter('datalist', 'datalist_Real.csv');   % e.g. 'datalist_Sham.csv' (stage F)           % start from the 2023 <id>_e.set (same trials as 2023) instead of the .cnt
p.parse(varargin{:});
o = p.Results;
src = data_dir;
setupPaths();

L = readtable(fullfile(o.outDir, o.datalist), 'TextType', 'string', 'Delimiter', ',');
L.evCodes = string(L.evCodes);
if ischar(which) || isstring(which)
    switch char(which)
        case 'pilot', L = L(strlength(L.pilot) > 0 & L.exclude == 0, :);
        case 'all',   L = L(L.exclude == 0, :);
        otherwise,    L = L(L.id == string(which), :);
    end
else
    L = L(ismember(L.id, string(which)), :);
end

ovFile = fullfile(o.outDir, 'ica_overrides.csv');
OV = table(zeros(0, 1), strings(0, 1), zeros(0, 1), zeros(0, 1), strings(0, 1), 'VariableNames', {'sub', 'cond', 'round', 'component', 'action'});
if exist(ovFile, 'file')
    R = readtable(ovFile, 'TextType', 'string', 'Delimiter', ',', 'CommentStyle', '#');
    if height(R) > 0
        OV = table(double(R.sub), strtrim(string(R.cond)), double(R.round), double(R.component), lower(strtrim(string(R.action))), ...
            'VariableNames', OV.Properties.VariableNames);
    end
end

lfFile = fullfile(o.outDir, 'label_overrides.csv');      % id,pulseSize,intensity,reason (EMG/VETA cross-check)
LF = table(strings(0, 1), zeros(0, 1), zeros(0, 1), 'VariableNames', {'id', 'pulseSize', 'intensity'});
if exist(lfFile, 'file'), LF = readtable(lfFile, 'TextType', 'string', 'Delimiter', ','); end
todo = false(height(L), 1); ovr = cell(height(L), 1);
for i = 1:height(L)
    ovr{i} = OV(OV.sub == L.sub(i) & OV.cond == L.cond(i), :);
    logF = fullfile(o.outDir, 'Logs', L.id(i) + "_procLog.mat");
    finF = fullfile(o.outDir, 'Final', L.id(i) + "_F_TEP.set");
    if o.force || ~exist(logF, 'file') || ~exist(finF, 'file')
        todo(i) = true;
    else
        S = load(logF); prev = S.procLog.overridesApplied;
        if isempty(prev), prev = OV([], :); end
        todo(i) = ~isequal(sortrows(prev), sortrows(ovr{i}));
    end
end
fprintf('%d of %d dataset(s) to run\n', sum(todo), height(L));
idx = find(todo)';
status = repmat(struct('id', "", 'ok', false, 'message', "", 'minutes', NaN), 1, numel(idx));
t23 = cell(height(L), 1);                                % trial2023 per paired pulse (checks/check_trials_vs_2023.m)
dtF = {fullfile(o.outDir, 'trials_vs_2023_detail.csv'), fullfile(o.outDir, 'trials_vs_2023_detail_Sham.csv')};
dtF = dtF(cellfun(@(f) exist(f, 'file') == 2, dtF));
if ~isempty(dtF)
    DT = [];
    for k = 1:numel(dtF), DT = [DT; readtable(dtF{k}, 'TextType', 'string')]; end %#ok<AGROW>
    for i = 1:height(L)
        d = DT(DT.id == L.id(i) & ~isnan(DT.newTrial), :); d = sortrows(d, 'newTrial'); t23{i} = d.trial2023';
    end
end
lfix = arrayfun(@(i) [LF.pulseSize(LF.id == L.id(i)) LF.intensity(LF.id == L.id(i))], (1:height(L))', 'uni', 0);
runOne = @(i) runDataset(L(i, :), ovr{i}, src, o.outDir, [o.extra, {'labelFix', lfix{i}, 'trial2023', t23{i}}], o.fromE);
if o.workers > 0
    if isempty(gcp('nocreate')), parpool('local', o.workers); end
    parfor k = 1:numel(idx), status(k) = runOne(idx(k)); end %#ok<PFBNS>
else
    for k = 1:numel(idx), status(k) = runOne(idx(k)); end
end
T = struct2table(status, 'AsArray', true);
writetable(T, fullfile(o.outDir, 'Logs', sprintf('batch_status_%s.csv', char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')))));
disp(T);
end

function s = runDataset(row, ovr, src, outDir, extra, fromE)
data_dir = 'path/to/data';   % release: set this path
setupPaths();
s = struct('id', row.id, 'ok', false, 'message', "", 'minutes', NaN);
t0 = tic;
try
    evCodes = num2cell(char(row.evCodes));
    evCodes = cellfun(@(c) c, evCodes, 'uni', 0);
    triallist = [];
    if row.labelSource == "triallist"
        tlDirs = {[data_dir '\spTMS_Formal_EMG, MEP\Real_TMS_VETA'], src};   % user 2026-10-10: triallists from G: first
        tlFile = '';
        for k = 1:numel(tlDirs), f = fullfile(tlDirs{k}, char(row.triallistFile)); if exist(f, 'file'), tlFile = f; break, end, end
        T = load(tlFile); triallist = T.triallist;
        fprintf('%s: triallist %s\n', row.id, tlFile);
    end
    maxNumCompThreads(3);
    preprocess_tep(char(row.id), evCodes, triallist, outDir, 'cntFile', fullfile(src, char(row.cnt)), ...
        'labelSource', char(row.labelSource), 'rate', row.rate, 'overrides', ovr, extra{:}, ...
        'eFile', fastifE(fromE, fullfile(src, [char(row.id) '_e.set'])));
    s.ok = true;
catch ME
    s.message = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
    fprintf(2, '%s FAILED: %s\n', row.id, ME.message);
end
s.minutes = toc(t0) / 60;
close all force
end

function setupPaths()
toolbox_dir = 'path/to/toolboxes';   % release: set this path
if isempty(which('pop_tesa_fastica'))
    addpath([toolbox_dir '\eeglab2021.1']);
    eeglab('nogui');
end
if isempty(which('fastica')), addpath([toolbox_dir '\eeglab2020_0\plugins\FastICA_25']); end
here = fileparts(mfilename('fullpath'));
addpath(here, '-begin');
end

function v = fastifE(c, f)
if c, v = f; else, v = ''; end
end
