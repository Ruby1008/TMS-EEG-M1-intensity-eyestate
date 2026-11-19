# TMS-EEG-M1-intensity-eyestate

Analysis code for the manuscript

**Stimulation intensity and eye state in M1 TMS-EEG: A sham-controlled study of cortical and corticospinal responses**

Li-An Chou, Yi-Chun Tsai, Hui-Yan Chiau, Wei-Kuang Liang, Chi-Hung Juan

Institute of Cognitive Neuroscience and Cognitive Intelligence and Precision Healthcare Research Center, National Central University, Taiwan

## The study

Real and sham single TMS pulses were delivered to the left primary motor cortex (M1) at 80%, 100% and 120% of the resting motor threshold (RMT), in separate eyes-open (EO) and eyes-closed (EC) blocks, while EEG and EMG were recorded. The code in this repository does the following:

- extracts motor-evoked potentials (MEPs) from the EMG of 23 participants;
- preprocesses the TMS-EEG of the 17 participants in the EEG analyses with EEGLAB and TESA (two FastICA rounds);
- decomposes the EEG with CEEMDAN. The TMS-evoked response is the empirical-mode reconstruction (ERM) of IMFs 4–7, and inter-trial phase coherence (ITPC) is computed per IMF;
- compares real and sham stimulation (GMFP against a bootstrap null);
- tests eye-state effects (pre-stimulus alpha power, surface Laplacian, cluster permutation test, TOST equivalence);
- fits single-trial mixed models that relate the P30 and pre-stimulus mu power to MEP amplitude;
- draws the figures of the main text and the supplementary material.

## Folder structure

```
TMS-EEG-M1-intensity-eyestate/
  preprocessing/            EEG preprocessing (MATLAB, EEGLAB + TESA, FastICA), AEP removal, QC table
  decomposition_ERM_ITPC/   CEEMDAN first layer, ERM (IMFs 4-7), ITPC (MATLAB and Python)
  EMG_MEP/                  MEP extraction from VETA trial tables, MEP rules, visual QC pages, shell drivers
  statistics/               TEP / real-vs-sham / eye-state / single-trial statistics (Python and R)
  figures/                  Figures 1-7, graphical abstract, supplementary figures, shared figure helper
  README.md  LICENSE  .gitignore
```

This repository holds code only. It contains no data, intermediate files or figures.

### Paths

Every script has one path variable near the top, or at the start of each MATLAB function that needs it:

- `data_dir` (`data_dir = 'path/to/data'` / `data_dir <- "path/to/data"` / `data_dir="path/to/data"`). Set it to the folder that holds the study data in the original project layout. The scripts expect the sub-folders `Revision_2026/...` (preprocessed data, intermediate results and outputs) and, for the raw EMG trial tables and the 2023 files, `spTMS_Formal_EMG, MEP/` and `spTMS_Formal_Sensor/` directly under `data_dir`.
- `toolbox_dir`, in the MATLAB files that load toolboxes. Set it to the folder that holds `eeglab2021.1`, `eeglab2020_0/plugins/FastICA_25` and `holo_move12_tb_20240530F` (see below).

In this repository the scripts are grouped by topic. Several of them import each other from their original folders (`sys.path.insert(...)`) or read each other's outputs at fixed relative locations. To re-run them, put each file back in its original folder under `data_dir`, as listed in the table at the end of this README.

## Software and versions

Versions were read from the analysis computer.

**MATLAB R2025a**, with:
- EEGLAB 2021.1, with the plugins TESA 1.1.1, ICLabel 1.4 and ANTeepimport 1.14 (reading the ANT `.cnt` files)
- FastICA 2.5 (`eeglab2020_0/plugins/FastICA_25`)
- `holo_move12_tb_20240530F`, a non-public toolbox that is not included here. It provides the CEEMDAN / Holo-Hilbert routines (e.g. `eemd1_tfpar`, `rcada_eiemdFD`, `mdtfx_aEEMD_D`) and the bundled SPM12 (version 6685) with FieldTrip. The FieldTrip version string is not recorded in the bundled copy; FieldTrip is used for the cluster permutation test in `erm_cluster_perm.m`.
- VETA (Jackson & Greenhouse), used upstream to convert the `.acq` EMG recordings into the VETA trial tables read by `mep_extract.m`. It is not included.

The non-public toolboxes (holo_move12_tb with CEEMDAN, SPM12 and FieldTrip) are named here only and are not redistributed.

**R 4.5.3** (2026-03-11, ucrt), with lme4 1.1.37, lmerTest 3.1.3, car 3.1.3, dplyr 1.1.4, ggplot2 4.0.2, cowplot 1.1.3, reshape2 1.4.4, extrafont 0.19 and R.matlab 3.7.0. The scripts and shell drivers call `Rscript`, so it must be on the PATH.

**Python 3.12.1**, with numpy 2.4.6, pandas 2.2.3, scipy 1.18.0, matplotlib 3.11.0, statsmodels 0.14.6, pingouin 0.6.1, openpyxl 3.1.5 and Pillow 12.2.0.

The TIFF conversion (`tiff_IJP.py`, `tiff_IJP_v3.py`) needs the sRGB and RSWOP.icm colour profiles from the standard Windows colour-profile folder `C:\Windows\System32\spool\drivers\color`.

The shell drivers (`EMG_MEP/*.sh`) were run with Git Bash on Windows.

## Order of the scripts

Supplementary items are numbered as in the submitted supplementary material (Figures S1–S10, Tables S1–S5). They were matched to the scripts by content.

### 1. EEG preprocessing (`preprocessing/`)
1. `check_trials_vs_2023.m`: trial key (trial2023) used to match EEG and MEP trials.
2. `run_preproc_batch.m`: batch driver. It calls `preprocess_tep.m`, which uses `pair_tms_intensity.m`, `decay_component_rule.m`, `icaOverviewFigure.m` and `tep_subject_qc.m`, for:
   - TESA pulse detection on FC3;
   - intensity labels and epoching;
   - removal of −5 to 10 ms;
   - PCA plus FastICA round 1 with the decay rule;
   - 0.1–80 Hz and 59–61 Hz filters;
   - FastICA round 2;
   - average reference.

   Outputs: `<id>_F_ERM.set` and `<id>_F_TEP.set`.
3. `export_preproc_summary.m`, then `build_qc_table.py` (data quality per condition, Table S1, 17-participant EEG sample), and `sa_tables_flowchart.py` (ICA statistics for Table S1).
4. `codes_without_pulse2.m`: list of trials without a detected pulse, used by the MEP rules.
5. `extract_prestim_phaseA_def.m`: pre-stimulus definition used by `statistics/erm_vs_tep_prestim.py`.
6. `aep_ica_review.m`, then `aep_remove.m`: ICA-based removal of auditory-evoked potentials (Figure S8).

### 2. CEEMDAN, ERM and ITPC (`decomposition_ERM_ITPC/`)
1. `hhsa_setup.m`: CEEMDAN parameters (50 noise pairs, noise level 0.2, 10 sifts, 11 IMFs, fixed random seed).
2. `run_layer1_batch.m` (calls `hhsa_layer1.m`), and `run_layer1_aepr.m` for the AEP-removed data.
3. `extract_all_real.m` (calls `erm_extract.m` and `itpc_extract.m`). The same two functions were run on each sham dataset from the command line; there is no saved driver for the sham runs.
4. `erm_single_trial_extract.m`: single-trial ERM P30 and N100.
5. `itpc_roi_export.m` (real) and `itpc_roi_export_1110.m` (sham).
6. `itpc_v2.py`: ITPC repeated-measures ANOVA, post hoc tests and FDR (Figure 5; ITPC rows of Tables S2 and S3). It loads the function definitions of `itpc_round15.py`. It is followed by `itpc_r15_prestim_check_v2.py`.
7. `ITPC_baseline_sens_IMF3to7.m`: ITPC baseline sensitivity analysis.
8. `n2_sham_itpc_1110.py`: ITPC with real and sham stimulation (Figure S9).
9. `imf_frequency_distribution.m`, then `fig4_bandwidth.py`: instantaneous-frequency distributions of the IMFs (Figure S3).

### 3. EMG / MEP (`EMG_MEP/`)
1. `mep_extract.m`: reads the VETA trial tables (`TMS_MEP_*_preprocessed.mat`).
2. `eeg_mep_trial_map.py`, `a1_022_mapping_check.py`, `c2_label_mapping_check.py` and `c2_cnt_compare.py`: EEG–MEP trial matching and intensity-label corrections. `c2_cnt_compare.py` reads the output of `c2_cnt_codes.m`, which is not included (see "Not included").
3. `make_qc_html.py` / `make_qc_html_v2.py`: visual MEP QC pages, then `qc_to_log_v2.py`: QC log.
4. `run_F1.sh` (MEP set `main_qc2`) and `run_F8.sh` / `run_F8b.sh` (MEP set `main_qc2_008`). These drivers run `mep_reprocess.py` with the fixed rules (RMS window −54 to −5 ms, no-pulse list, intensity overrides, QC log) and several downstream scripts. `mep_reprocess.py` also writes the MEP probability matrix used by the graphical abstract. Some steps in the drivers call scripts that are not in this release (see "Not included"); those steps can be skipped.
5. `mep_stats_A_vs_B_reprocessed.py`: MEP repeated-measures ANOVAs on log and raw amplitude, with normality and non-parametric checks (Results, Figure 2 brackets). `validate_anova_vs_pingouin.py` checks the ANOVA function against pingouin.

### 4. Statistics (`statistics/`)

**Conventional TEP:**
- `tep_roi_waveforms.m` and `tep_allch_extract.py`, then `tep_windows.py`, then `tep_stats.py` and `tep_supp_measures.py` (sliding-window and component ANOVAs, band-pass control, GMFP).
- `early_component_dependence.py`.
- Shared helpers: `tep_io.py` and `tep_stats_lib.py`.

**Real vs sham:**
- `r21_sham.py` (GMFP 15–60 ms against a 2000-draw bootstrap null), then `r21_round15.py` (local-baseline null; Figure 4).
- Then `R2_r21_round15_noOverride_rev1108.py`, `n3_reliable_windows_1110.py` (scalp distribution at 80% RMT, Figure S5) and `r21_n100_confound.py` (N100, real − sham).

**Pre-stimulus interval:** `erm_vs_tep_prestim.py` and `erm_prestim_variants.py` (pre-stimulus ERM and band-pass control, Figure S4).

**AEP control:** `aep_compare_tep.py`, `aep_components.py` and `aep_compare_erm_itpc.py` (Figure S8).

**Eye state:**
- `prestim_alpha_check.py`, then `round9_extract.py` (alpha power, surface Laplacian; Figure 6), then `eyestate_null_round8.py` (calls `eyestate_null_round8.R`) and `eyestate_itpc_v2.py` (calls `eyestate_itpc_v2.R`): eye-state effects, Table S2.
- `eyestate_mep_rows_reprocessed.py`: Table S2 MEP rows. It writes and runs a temporary R script.
- `erm_cluster_perm.m`: FieldTrip cluster permutation test, 5000 permutations (Table S2).
- `tost_eo_ec_reprocessed.py`: equivalence of EO and EC, paired TOST (Table S3, Figure 6).
- `n1_conv_tep_eye_1110.py`, `n4_residual_components_1110.py` and `n4_lmm_1110.R`. Shared helper: `r1110_lib.py`.

**Single trial (P30, mu power, MEP):**
- `single_trial_robustness_extract.py`, then `build_single_trial_inputs_reprocessed.py`, then `single_trial_robustness_reprocessed.R` (Figure 7, Figure S10).
- Then `TableS5_fullmodel.R`: single-trial models of log MEP on the P30 split into within- and between-participant components, without and with one additional covariate. These are the models of Table S4; the script is named after the earlier table number.
- Then `TableS5_randomslope.R`: random-slope versions of the Table S4 and S5 models. It reads the output of `TableS5_fullmodel.R`.
- Mu power and P30:
  - `mu_3e_reprocessed.R`: main rows of Table S5;
  - `n5_extract_1110.py`, then `n5_lmm_1110.R`, then `n5_maps_1110.py`: alternative definitions in Table S5;
  - `round9_mu_stats.py`: channel-wise models;
  - `mu_P30_partE.R`: per-participant slopes, Figure 7;
  - `mu_quintiles_reprocessed.py`: Figure 7.

**Checkpoints (re-run on 2026-11-18, MEP set `main_qc2_008`, 17 participants):**
- `TableS5_fullmodel.R`, Table S4 main models:
  - 100% RMT: P30_within b = 0.0327 (SE 0.0069), p < .001; P30_between b = 0.0189 (0.1437), p = .897.
  - 120% RMT: P30_within b = 0.0115 (0.0045), p = .011; P30_between b = 0.1514 (0.1056), p = .170.
- `mu_3e_reprocessed.R`, Table S5 mu models of log MEP:
  - 100% RMT: b = −0.101 (0.072), p = .165.
  - 120% RMT: b = −0.028 (0.067), p = .675.

### 5. Figures (`figures/`)

| Figure | Script |
|---|---|
| Figure 1 | `fig1_1114.py` |
| Figure 2 | `Fig2_1114.R` (MEPSET=main_qc2) |
| Figure 3 | `fig3_1114.py` |
| Figure 4 | `fig4_IJP.py`, TIFF export `tiff_IJP.py` |
| Figure 5 | `fig5_1114.py` |
| Figure 6 | `fig6_1114.py` (MEPSET=main_qc2) |
| Figure 7 | `fig7_IJP_v3.py` (MEPSET=main_qc2_008), TIFF export `tiff_IJP_v3.py` |
| Graphical abstract | `ga_IJP_v2.py` |
| Figure S1 Preprocessing pipeline | `figS2_1109i.py` |
| Figure S2 Single-participant TEPs | `figS3_1108t.py`; PDF of single-participant TEPs: `sd_single_subject_tep_final_b_potential.py` |
| Figure S3 Instantaneous-frequency distributions of the IMFs | `figS5_oldS4_1108n_v2.m` (with `decomposition_ERM_ITPC/imf_frequency_distribution.m`) |
| Figure S4 Pre-stimulus interval of the ERM | `figS4_oldS7_1108n.py` |
| Figure S5 Scalp distribution of the early response at 80% RMT | `statistics/n3_reliable_windows_1110.py` |
| Figure S6 ERM at the ROI, real vs sham; Figure S7 ERM at Cz, real vs sham | `figS6_S8_oldS15_S14_1108n.py` |
| Figure S8 ERM before and after AEP removal | `figS12_oldS21_1108n.py` |
| Figure S9 ITPC with real and sham stimulation | `decomposition_ERM_ITPC/n2_sham_itpc_1110.py` |
| Figure S10 Per-participant relation between single-trial P30 and log MEP | `figS10_1108v.py` (MEPSET=main_qc2_008) |
| Trial counts (Table S1) | `trial_counts_1108n.py` |

`rev1108_helpers.py` is a shared helper imported by several of the supplementary figure scripts.

The numbers in the file names (e.g. `figS2_`, `oldS7`, `S14`) are earlier figure numbers. The file names were kept unchanged.

## Not included

- **Excluded for privacy:** `c2_cnt_codes.m` and `c2_c1_summary.py`. These two EEG–MEP trial-matching helpers contain participant name initials in raw file names, so they are not included in any form.
- **Items that are not in the submitted supplementary material** (they existed only in an earlier version). Their scripts are left out because their numbers could not be checked against the submitted text:
  - MEP probability table: `make_table_mep_probability.py`;
  - grand-average conventional TEPs at 61 channels: `figS1_1108n.py`;
  - pre-stimulus power spectra: `figS11_1108t.py`;
  - late intensity effects: `late_window_imf_cz.py`, `late_window_mep_covariate_reprocessed.R`;
  - cross-participant TEP–MEP correlation table: `fig6_input_reprocessed.py`, `Fig6_new_TEP_MEP.R`, `fig6_pearson_logmep_reprocessed.py`.
- **Not part of this release:**
  - `joint_lmm_1113.R`, a second implementation of the Table S4 joint-covariate model. It gives the same estimates as `TableS5_fullmodel.R` but needs outputs of an unfinalized analysis.
  - `fig7C_randomslope.R` and `fig7C_export_trials.py`: random-slope channel map of an earlier version of Figure 7, not in the final figure.
- **Analyses not finalized for the manuscript, not included:**
  - raw vs log MEP comparison (`logMEP_2026-11-12`);
  - pre-stimulus Figure 3A variants (`Fig3A_prestim_2026-11-14`);
  - `figS3_1110.py`, `itpc_fig78_G.py`, `itpc_posthoc_pooled_prestim.py`, `R1_itpc_v2_tables_rev1108.py`, `g_checks_1110.py`, `r1_1_main17.py`, `tep_supp_lmm.R`;
  - the GG-df export scripts.
- **Called by the shell drivers but not included:** `Fig2_MEP_boxplot_a3.R`, `Fig2_MEP_boxplot_reprocessed.R`, `Fig2_p6_qc2.R`, `Fig6_new_TEP_MEP.R`, `fig6_input_reprocessed.py`, `fig6_pearson_logmep_reprocessed.py`, `make_table_mep_probability.py`, `late_window_mep_covariate_reprocessed.R`, `bootstrap_itpc_mep_17_reprocessed.R`, `fig6_single_trial_combined_reprocessed.R`, `fig6_single_trial_lmm_reprocessed.R`, `fig6_single_trial_supp_reprocessed.R`, `prestim_power_lmm_reprocessed.R`, `prestim_power_lmm_v2_reprocessed.R`, `corr_spec_ladder_reprocessed.py`, `fig6_a3.py`, `fig6_reprocessed.py`, `fig7_a3.py`, `fig7_reprocessed.py`, `single_trial_per_participant_plot_reprocessed.py`, `tep_mep_sliding_reprocessed.py`. These are earlier figure versions or analyses that are not in the submitted manuscript.
- **Optional and not used in the main run:** `epochs_from_e.m`, called by `preprocess_tep.m` only with the option `fromE`.
- **Other files:** electrode positions (`waveguard64_woRef.ced`, read from `data_dir`) are a data file and are not included.

## Code development

The scripts in this repository were written with the assistance of an AI coding assistant (Claude, Anthropic) and reviewed by the authors.

## Data availability

The data are not publicly available due to ethical restrictions on participant privacy. De-identified data are available from the corresponding author (Chi-Hung Juan, chijuan@cc.ncu.edu.tw) on reasonable request, subject to approval by the Ethics Committee of National Taiwan University (201910EM003).

## Citation

Citation: to be added after publication.

## License

MIT License (see `LICENSE`).

## Original location of each file

Paths are relative to `data_dir`.

| File in this repository | Original folder |
|---|---|
| `preprocessing/preprocess_tep.m` | `Revision_2026/Preproc/code/` |
| `preprocessing/run_preproc_batch.m` | `Revision_2026/Preproc/code/` |
| `preprocessing/pair_tms_intensity.m` | `Revision_2026/Preproc/code/` |
| `preprocessing/decay_component_rule.m` | `Revision_2026/Preproc/code/` |
| `preprocessing/icaOverviewFigure.m` | `Revision_2026/Preproc/code/` |
| `preprocessing/tep_subject_qc.m` | `Revision_2026/Preproc/code/` |
| `preprocessing/export_preproc_summary.m` | `Revision_2026/Preproc/code/` |
| `preprocessing/build_qc_table.py` | `Revision_2026/Preproc/code/` |
| `preprocessing/extract_prestim_phaseA_def.m` | `Revision_2026/Preproc/code/` |
| `preprocessing/check_trials_vs_2023.m` | `Revision_2026/Preproc/code/checks/` |
| `preprocessing/codes_without_pulse2.m` | `Revision_2026/Manuscript_Figures/code/` |
| `preprocessing/aep_ica_review.m` | `Revision_2026/Manuscript_Figures/code/` |
| `preprocessing/aep_remove.m` | `Revision_2026/Manuscript_Figures/code/` |
| `preprocessing/sa_tables_flowchart.py` | `Revision_2026/Manuscript_Figures/code/` |
| `decomposition_ERM_ITPC/hhsa_setup.m` | `Revision_2026/Manuscript_Figures/code/` |
| `decomposition_ERM_ITPC/hhsa_layer1.m` | `Revision_2026/Manuscript_Figures/code/` |
| `decomposition_ERM_ITPC/run_layer1_batch.m` | `Revision_2026/Manuscript_Figures/code/` |
| `decomposition_ERM_ITPC/run_layer1_aepr.m` | `Revision_2026/Manuscript_Figures/code/` |
| `decomposition_ERM_ITPC/erm_extract.m` | `Revision_2026/Manuscript_Figures/code/` |
| `decomposition_ERM_ITPC/itpc_extract.m` | `Revision_2026/Manuscript_Figures/code/` |
| `decomposition_ERM_ITPC/extract_all_real.m` | `Revision_2026/Manuscript_Figures/code/` |
| `decomposition_ERM_ITPC/erm_single_trial_extract.m` | `Revision_2026/Manuscript_Figures/code/` |
| `decomposition_ERM_ITPC/itpc_roi_export.m` | `Revision_2026/Manuscript_Figures/code/` |
| `decomposition_ERM_ITPC/itpc_v2.py` | `Revision_2026/Manuscript_Figures/code/` |
| `decomposition_ERM_ITPC/itpc_round15.py` | `Revision_2026/Manuscript_Figures/code/` |
| `decomposition_ERM_ITPC/itpc_r15_prestim_check_v2.py` | `Revision_2026/Manuscript_Figures/code/` |
| `decomposition_ERM_ITPC/ITPC_baseline_sens_IMF3to7.m` | `Revision_2026/Manuscript_Figures/code/` |
| `decomposition_ERM_ITPC/imf_frequency_distribution.m` | `Revision_2026/Manuscript_Figures/code/` |
| `decomposition_ERM_ITPC/fig4_bandwidth.py` | `Revision_2026/Manuscript_Figures/code/` |
| `decomposition_ERM_ITPC/itpc_roi_export_1110.m` | `Revision_2026/Manuscript_Figures/Review_2026-11-10/code/` |
| `decomposition_ERM_ITPC/n2_sham_itpc_1110.py` | `Revision_2026/Manuscript_Figures/Review_2026-11-10/code/` |
| `EMG_MEP/mep_extract.m` | `Revision_2026/Manuscript_Figures/code/mep_reprocess/` |
| `EMG_MEP/mep_reprocess.py` | `Revision_2026/Manuscript_Figures/code/mep_reprocess/` |
| `EMG_MEP/validate_anova_vs_pingouin.py` | `Revision_2026/Manuscript_Figures/code/mep_reprocess/` |
| `EMG_MEP/eeg_mep_trial_map.py` | `Revision_2026/Manuscript_Figures/code/` |
| `EMG_MEP/make_qc_html.py` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/` |
| `EMG_MEP/a1_022_mapping_check.py` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/` |
| `EMG_MEP/c2_label_mapping_check.py` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/` |
| `EMG_MEP/c2_cnt_compare.py` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/` |
| `EMG_MEP/run_pipeline.sh` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/` |
| `EMG_MEP/make_qc_html_v2.py` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/v2/` |
| `EMG_MEP/qc_to_log_v2.py` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/v2/` |
| `EMG_MEP/run_F1.sh` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/v2/` |
| `EMG_MEP/run_F8.sh` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/v2/` |
| `EMG_MEP/run_F8b.sh` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/v2/` |
| `EMG_MEP/mep_stats_A_vs_B_reprocessed.py` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/v2/code008/` |
| `statistics/r21_sham.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/r21_round15.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/r21_n100_confound.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/tep_roi_waveforms.m` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/tep_allch_extract.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/tep_windows.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/tep_stats.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/tep_supp_measures.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/tep_io.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/tep_stats_lib.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/early_component_dependence.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/erm_vs_tep_prestim.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/erm_prestim_variants.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/aep_compare_tep.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/aep_components.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/aep_compare_erm_itpc.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/prestim_alpha_check.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/round9_extract.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/round9_mu_stats.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/eyestate_null_round8.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/eyestate_null_round8.R` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/eyestate_itpc_v2.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/eyestate_itpc_v2.R` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/erm_cluster_perm.m` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/single_trial_robustness_extract.py` | `Revision_2026/Manuscript_Figures/code/` |
| `statistics/R2_r21_round15_noOverride_rev1108.py` | `Revision_2026/Manuscript_Figures/Figures_all_2026-11-08/code/` |
| `statistics/n3_reliable_windows_1110.py` | `Revision_2026/Manuscript_Figures/Review_2026-11-10/code/` |
| `statistics/n1_conv_tep_eye_1110.py` | `Revision_2026/Manuscript_Figures/Review_2026-11-10/code/` |
| `statistics/n4_residual_components_1110.py` | `Revision_2026/Manuscript_Figures/Review_2026-11-10/code/` |
| `statistics/n4_lmm_1110.R` | `Revision_2026/Manuscript_Figures/Review_2026-11-10/code/` |
| `statistics/r1110_lib.py` | `Revision_2026/Manuscript_Figures/Review_2026-11-10/code/` |
| `statistics/n5_extract_1110.py` | `Revision_2026/Manuscript_Figures/Review_2026-11-10/code/` |
| `statistics/n5_lmm_1110.R` | `Revision_2026/Manuscript_Figures/Review_2026-11-10/code/` |
| `statistics/n5_maps_1110.py` | `Revision_2026/Manuscript_Figures/Review_2026-11-10/code/` |
| `statistics/tost_eo_ec_reprocessed.py` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/v2/code008/` |
| `statistics/eyestate_mep_rows_reprocessed.py` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/v2/code008/` |
| `statistics/build_single_trial_inputs_reprocessed.py` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/v2/code008/` |
| `statistics/single_trial_robustness_reprocessed.R` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/v2/code008/` |
| `statistics/TableS5_fullmodel.R` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/v2/code008/` |
| `statistics/TableS5_randomslope.R` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/v2/code008/` |
| `statistics/mu_3e_reprocessed.R` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/v2/code008/` |
| `statistics/mu_P30_partE.R` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/v2/code008/` |
| `statistics/mu_quintiles_reprocessed.py` | `Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/v2/code008/` |
| `figures/fig1_1114.py` | `Revision_2026/Manuscript_Figures/Figures_2026-11-14_minor/code/` |
| `figures/Fig2_1114.R` | `Revision_2026/Manuscript_Figures/Figures_2026-11-14_minor/code/` |
| `figures/fig3_1114.py` | `Revision_2026/Manuscript_Figures/Figures_2026-11-14_minor/code/` |
| `figures/fig4_IJP.py` | `Revision_2026/Manuscript_Figures/Figure4_2026-11-11/code/` |
| `figures/tiff_IJP.py` | `Revision_2026/Manuscript_Figures/Figures_IJP_2026-11-10/code/` |
| `figures/fig5_1114.py` | `Revision_2026/Manuscript_Figures/Figures_2026-11-14_minor/code/` |
| `figures/fig6_1114.py` | `Revision_2026/Manuscript_Figures/Figures_2026-11-14_minor/code/` |
| `figures/fig7_IJP_v3.py` | `Revision_2026/Manuscript_Figures/Figures_2026-11-14_minor/code/` |
| `figures/tiff_IJP_v3.py` | `Revision_2026/Manuscript_Figures/Figures_IJP_2026-11-10/code/` |
| `figures/ga_IJP_v2.py` | `Revision_2026/Manuscript_Figures/Figures_IJP_2026-11-10/code/` |
| `figures/figS2_1109i.py` | `Revision_2026/Manuscript_Figures/Figures_2026-11-09i/code/` |
| `figures/figS3_1108t.py` | `Revision_2026/Manuscript_Figures/Figures_2026-11-08t/code/` |
| `figures/sd_single_subject_tep_final_b_potential.py` | `Revision_2026/Manuscript_Figures/FigS_SD_2026-11-08/final_b_potential/code/` |
| `figures/figS5_oldS4_1108n_v2.m` | `Revision_2026/Manuscript_Figures/Figures_2026-11-08n/code/` |
| `figures/figS4_oldS7_1108n.py` | `Revision_2026/Manuscript_Figures/Figures_2026-11-08n/code/` |
| `figures/figS6_S8_oldS15_S14_1108n.py` | `Revision_2026/Manuscript_Figures/Figures_2026-11-08n/code/` |
| `figures/figS12_oldS21_1108n.py` | `Revision_2026/Manuscript_Figures/Figures_2026-11-08n/code/` |
| `figures/figS10_1108v.py` | `Revision_2026/Manuscript_Figures/Figures_2026-11-08v/code/` |
| `figures/trial_counts_1108n.py` | `Revision_2026/Manuscript_Figures/Figures_2026-11-08n/code/` |
| `figures/rev1108_helpers.py` | `Revision_2026/Manuscript_Figures/Figures_all_2026-11-08/code/` |
