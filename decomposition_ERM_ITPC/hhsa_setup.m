function P = hhsa_setup()
% HHSA_SETUP  Paths and fixed parameters for re-running the manuscript's EMD/HHSA pipeline
% on the 2026 preprocessed data (Revision_2026\Preproc\Final\*_F_ERM.set).
% Parameters are copied from the D.history of the original result files in
% <data_dir>/spTMS_Formal_Sensor (e.g. timf_spm_spTMS_001_Real_EO_F.mat); the only addition is a
% fixed random seed (the original runs used an unseeded randn).
data_dir = 'path/to/data';   % release: set this path
toolbox_dir = 'path/to/toolboxes';   % release: set this path
P.tb      = [toolbox_dir '\holo_move12_tb_20240530F\holo_move12_tb'];
P.eeglab  = [toolbox_dir '\eeglab2021.1'];
P.src     = [data_dir '\Revision_2026\Preproc\Final'];
P.out     = [data_dir '\Revision_2026\Manuscript_Figures\HHSA'];
P.suffix  = '_F_ERM';                        % input file = <id><suffix>.set in P.src
% AEP-ICA robustness (2026-10-15): setenv('HHSA_VARIANT','AEPr') reads Preproc_AEPr\Final\<id>_F_ERM_AEPr.set and
% writes to Manuscript_Figures\HHSA_AEPr (same parameters and seed); unset = main analysis (unchanged)
if strcmp(getenv('HHSA_VARIANT'), 'AEPr')
    P.src = [data_dir '\Revision_2026\Preproc_AEPr\Final'];
    P.out = [data_dir '\Revision_2026\Manuscript_Figures\HHSA_AEPr'];
    P.suffix = '_F_ERM_AEPr';
end
P.seed    = 20261012;                        % rng seed, set before every eemd1_tfpar / aeemd2layer call
P.dropChan = {'M1','M2'};                    % mastoids: not in the original 61-channel SPM files
% first layer (eemd1_tfpar) - original: eetype 2 (rcada_eiemdFD), NE 50, noise 0.2, shift 1, [-2000 1999] ms, quadm
P.l1 = struct('channels', {{'EEG'}}, 'DCMEMD', 0, 'eetype', 2, 'Flandrin', 0, 'parallel', 1, 'CoreN', 5, ...
    'timewin', [-2000 1999], 'shiftLevel', 1, 'NEnsemble', 50, 'ENoise', 0.2, 'NEnsemble2', 20, 'ENoise2', 0.2, ...
    'UpSampX2', 0, 'phase', 1, 'fmExp', 0, 'ifmethod', 'quadm', 'from1stEmd', 0, 'nyquist', 0);
% second layer (aeemd2layer_tfpar_SV) - original: from1stEmd 1, NE2 20, noise 0.2, sdownLV 2, shift 0
P.l2 = struct('channels', {{'EEG'}}, 'DCMEMD', 0, 'eetype', 2, 'Flandrin', 0, 'parallel', 1, 'CoreN', 5, ...
    'from1stEmd', 1, 'keep1stfm', 0, 'PartImfs', 0, 'phase', 1, 'phase2', 1, 'NEnsemble2', 20, 'ENoise2', 0.2, ...
    'NEnsemble', 20, 'ENoise', 0.2, 'ifmethod', 'quadm', 'ifmethod2', 'quadm', 'sdownLV', 2, 'timewin', [-2000 1999], 'shiftLevel', 0);
addpath(fullfile(P.tb, 'spm12')); addpath(P.tb);
for d = {'EEMD', 'emd', '2layer_EMD', 'conn', 'cubicSpline', 'statistics', 'RMANOVA'}, addpath(fullfile(P.tb, d{1})); end   % helpers used on the parfor workers (e.g. EEMD\FAqzc.m)
addpath(P.eeglab);
if ~exist(P.out, 'dir'), mkdir(P.out); end
end
