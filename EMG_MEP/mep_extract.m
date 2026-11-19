function mep_extract(in_dir, pattern, out_dir)
% MEP_EXTRACT  Export raw EMG segments around each TMS pulse from VETA files.
%
%   mep_extract(in_dir, pattern, out_dir)
%
%   Reads every VETA trials table matching fullfile(in_dir, pattern)
%   (e.g. '*_preprocessed.mat' or '*_acq_to_Veta*.mat'; the raw EMG in
%   trials.ch4 is identical in both) and writes two CSV files to out_dir:
%
%     mep_trials_meta.csv   one row per trial: file, subject, eye, trial,
%                           TMS_Condition, RMT, IPI, pulse sample, and the
%                           old VETA values (trial_accept, amplitude, RMS,
%                           artloc, MEP time) for comparison only
%     mep_waveforms.csv     one row per trial: EMG from -100 to +150 ms
%                           around the pulse (251 samples, 1000 Hz, mV)
%
%   The pulse is the first sample of each epoch at which the trigger
%   channel (trials.ch2) exceeds half of its range. No filtering, no
%   detrending and no event detection are done here; all rules are applied
%   in mep_reprocess.py.
%
%   Read-only with respect to the input files.

if nargin < 3
    error('usage: mep_extract(in_dir, pattern, out_dir)');
end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

fs = 1000;                          % Hz, checked against parameters.sampling_rate
pre = 100; post = 150;              % ms around the pulse
files = dir(fullfile(in_dir, pattern));
if isempty(files), error('no files match %s', fullfile(in_dir, pattern)); end

meta = {};
wave = zeros(0, pre + post + 1);
for f = 1:numel(files)
    fn = files(f).name;
    S = load(fullfile(files(f).folder, fn));
    if isfield(S, 'parameters') && isfield(S.parameters, 'sampling_rate') && S.parameters.sampling_rate ~= fs
        error('%s: sampling rate %g, expected %g', fn, S.parameters.sampling_rate, fs);
    end
    T = S.trials;
    tok = regexp(fn, 'TMS_MEP_(\d{3})_\d+_(EO|EC)', 'tokens', 'once');
    if isempty(tok), error('cannot parse subject/eye from %s', fn); end
    stim = 'Real';
    if contains(fn, 'Sham', 'IgnoreCase', true), stim = 'Sham'; end
    vars = T.Properties.VariableNames;
    for i = 1:height(T)
        emg = T.ch4{i}(:)';
        trg = T.ch2{i}(:)';
        thr = min(trg) + (max(trg) - min(trg)) / 2;
        k = find(trg > thr, 1);
        if isempty(k) || (max(trg) - min(trg)) < 1
            k = NaN;
            seg = nan(1, pre + post + 1);
        else
            idx = (k - pre):(k + post);
            seg = nan(1, numel(idx));
            ok = idx >= 1 & idx <= numel(emg);
            seg(ok) = emg(idx(ok));
        end
        meta(end + 1, :) = {fn, tok{1}, tok{2}, stim, getv(T, vars, 'trial_number', i), ...
            getv(T, vars, 'TMS_Condition', i), getv(T, vars, 'RMT', i), getv(T, vars, 'IPI', i), k, ...
            getv(T, vars, 'trial_accept', i), getv(T, vars, 'ch4_MEP_amplitude', i), ...
            getv(T, vars, 'ch4_RMS_preMEP', i), getv(T, vars, 'artloc', i), getv(T, vars, 'ch4_MEP_time', i)}; %#ok<AGROW>
        wave(end + 1, :) = seg; %#ok<AGROW>
    end
    fprintf('%s: %d trials\n', fn, height(T));
end

M = cell2table(meta, 'VariableNames', {'file', 'subject', 'eye', 'stim', 'trial_number', ...
    'TMS_Condition', 'RMT', 'IPI', 'pulse_sample', 'veta_trial_accept', 'veta_amplitude', ...
    'veta_rms_pre', 'veta_artloc', 'veta_mep_time'});
M.row_id = (1:height(M))';
writetable(M, fullfile(out_dir, 'mep_trials_meta.csv'));

W = array2table(wave, 'VariableNames', arrayfun(@(t) sprintf('t%+04d', t), -pre:post, 'UniformOutput', false));
W = [table((1:height(W))', 'VariableNames', {'row_id'}), W];
writetable(W, fullfile(out_dir, 'mep_waveforms.csv'));
fprintf('wrote %d trials from %d files to %s\n', height(M), numel(files), out_dir);
end

function v = getv(T, vars, name, i)
if any(strcmp(vars, name)), v = T.(name)(i); else, v = NaN; end
end
