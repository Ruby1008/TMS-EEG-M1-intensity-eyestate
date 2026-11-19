function run_layer1_aepr()
% RUN_LAYER1_AEPR  First-layer CEEMDAN on the AEP-removed real-stimulation data (AEP-ICA robustness, HANDOFF section 15),
% only the 17 main-analysis participants x EO/EC (34 datasets), then ERM and ITPC extraction (erm_extract, itpc_extract).
% Same parameters/seed as the main analysis (hhsa_setup); inputs Preproc_AEPr\Final\<id>_F_ERM_AEPr.set,
% outputs Manuscript_Figures\HHSA_AEPr\. Resumable (existing timf_ / erm_ / itpc_ files are skipped).
% Log: HHSA_AEPr\layer1_batch.log
setenv('HHSA_VARIANT', 'AEPr'); P = hhsa_setup();
subs = [1 2 3 7 8 9 12 15 16 17 19 21 22 24 26 27 28];      % 17 with complete real-stimulation EEG (P\QC\QC_summary_wide.csv)
ids = {}; for s = subs, for e = {'EO', 'EC'}, ids{end + 1} = sprintf('spTMS_%03d_Real_%s', s, e{1}); end, end %#ok<AGROW>
logF = fullfile(P.out, 'layer1_batch.log'); fid = fopen(logF, 'a'); fprintf(fid, '%s\tAEPr batch: %d datasets\n', datestr(now), numel(ids)); fclose(fid);
for k = 1:numel(ids)
    fid = fopen(logF, 'a');
    try
        info = hhsa_layer1(ids{k});
        if ~exist(fullfile(P.out, ['erm_' ids{k} '.mat']), 'file'), erm_extract(ids{k}); end
        if ~exist(fullfile(P.out, ['itpc_' ids{k} '.mat']), 'file'), itpc_extract(ids{k}); end
        if info.skipped, fprintf(fid, '%s\t%s\tlayer 1 skipped (exists); erm/itpc ok\n', datestr(now), ids{k});
        else, fprintf(fid, '%s\t%s\tchan %d\ttrials %d\t%.1f min\n', datestr(now), ids{k}, info.nChan, info.nTrials, info.minutes); end
    catch ME
        fprintf(fid, '%s\t%s\tFAILED: %s\n', datestr(now), ids{k}, strrep(getReport(ME, 'basic', 'hyperlinks', 'off'), newline, ' | '));
        for pre = {'timf_', 'timph_', 'timfm_', 'timam_'}
            for ext = {'.mat', '.dat'}, f = fullfile(P.out, [pre{1} 'spm_' ids{k} ext{1}]); if exist(f, 'file'), delete(f); end, end
        end
    end
    fclose(fid);
end
end
