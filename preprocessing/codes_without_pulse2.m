% Codes grouped into trials (codes are repeated within a trial): a new trial starts when the gap between
% consecutive codes > 2 s. For each code-trial: pulse within [-5 500] ms of its first code? -> trial index in
% stimulation order (= MEP trial number) of the trials without a pulse. Also pulses without any code-trial.
data_dir = 'path/to/data';   % release: set this path
toolbox_dir = 'path/to/toolboxes';   % release: set this path
addpath([toolbox_dir '\eeglab2021.1']); eeglab nogui;
src = data_dir; out = [data_dir '\Revision_2026\Manuscript_Figures\data'];
L = readtable([data_dir '\Revision_2026\Preproc\datalist_Real.csv'], 'TextType', 'string', 'Delimiter', ',');
ids = ["spTMS_012_Real_EO", "spTMS_012_Real_EC", "spTMS_015_Real_EC", "spTMS_021_Real_EO", "spTMS_022_Real_EO"]; rows = {};
for id = ids
    i = find(L.id == id); C = pop_loadeep_v4(fullfile(src, char(L.cnt(i))));
    C = pop_tesa_findpulse(C, 'FC3', 'refract', 3, 'rate', L.rate(i), 'tmsLabel', 'TMS', 'plots', 'off');
    ty = cellfun(@num2str, {C.event.type}, 'uni', 0); lat = [C.event.latency];
    cLat = lat(ismember(ty, {'1', '2', '3'})); cTy = ty(ismember(ty, {'1', '2', '3'})); [cLat, o] = sort(cLat); cTy = cTy(o);
    st = [1 find(diff(cLat) > 2000) + 1]; en = [st(2:end) - 1 numel(cLat)];
    pLat = sort(lat(strcmp(ty, 'TMS'))); pLat = pLat([true, diff(pLat) > 4000]);
    noP = []; usedP = false(size(pLat));
    for k = 1:numel(st)
        d = pLat - cLat(st(k)); j = find(d >= -5 & d <= 500, 1);
        if isempty(j), noP(end + 1) = k; else, usedP(j) = true; end %#ok<AGROW>
    end
    fprintf('%s: code-trials %d, pulses %d, trials without pulse %s (codes %s), pulses without code %d\n', id, numel(st), numel(pLat), mat2str(noP), strjoin(cTy(st(noP)), ','), sum(~usedP));
    for k = noP, rows(end + 1, :) = {id, k, cTy{st(k)}, numel(st), numel(pLat), round(cLat(st(k)))}; end %#ok<AGROW>
end
writetable(cell2table(rows, 'VariableNames', {'id', 'stimOrderIndex', 'code', 'nCodeTrials', 'nPulses', 'codeLatency_ms'}), fullfile(out, 'trials_without_pulse.csv'));
