function run_layer1_batch(which)
% RUN_LAYER1_BATCH  First-layer CEEMDAN (hhsa_layer1) over the processed datasets, resumable.
%   run_layer1_batch('Real')   all Real datasets with Final\<id>_F_ERM.set (stage-D excluded ones too: they
%                              are needed for the Supplementary); datasets with an existing timf_ file are skipped
%   run_layer1_batch('Sham')
% Log: HHSA\layer1_batch.log (one line per dataset). ~26 min per dataset with 5 workers.
if nargin < 1, which = 'Real'; end
P = hhsa_setup();
F = dir(fullfile(P.src, ['spTMS_*_' which '_*_F_ERM.set'])); ids = erase({F.name}, '_F_ERM.set');
logF = fullfile(P.out, 'layer1_batch.log'); fid = fopen(logF, 'a');
fprintf(fid, '%s\tbatch %s: %d datasets\n', datestr(now), which, numel(ids)); fclose(fid);
for k = 1:numel(ids)
    fid = fopen(logF, 'a');
    try
        info = hhsa_layer1(ids{k});
        if info.skipped, fprintf(fid, '%s\t%s\tskipped (exists)\n', datestr(now), ids{k});
        else, fprintf(fid, '%s\t%s\tchan %d\ttrials %d\t%.1f min\n', datestr(now), ids{k}, info.nChan, info.nTrials, info.minutes); end
    catch ME
        fprintf(fid, '%s\t%s\tFAILED: %s\n', datestr(now), ids{k}, strrep(getReport(ME, 'basic', 'hyperlinks', 'off'), newline, ' | '));
        % remove partial outputs so the dataset is retried on the next run
        for pre = {'timf_', 'timph_', 'timfm_', 'timam_'}
            for ext = {'.mat', '.dat'}, f = fullfile(P.out, [pre{1} 'spm_' ids{k} ext{1}]); if exist(f, 'file'), delete(f); end, end
        end
    end
    fclose(fid);
end
end
