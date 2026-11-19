#!/usr/bin/env bash
data_dir="path/to/data"   # release: set this path
# Part 8 item 1 (user 2026-10-23): 008 EC in the single-trial analyses. Same steps as qc/run_pipeline.sh (unchanged) with the step-7 programs taken
# from qc/v2/code008/ (make_code008.py: exclusion of spTMS_008_Real_EC removed). New output folders MEP_reprocess_2026-10-22/<SET>_008/; existing folders untouched.
set -u
R="${data_dir}\Revision_2026\Manuscript_Figures\MEP_reprocess_2026-10-22"
RU="${data_dir}/Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22"
RS="Rscript"
V="${data_dir}/Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/v2"
W="${data_dir}/Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc"
C8="$RU/qc/v2/code008"
export PYTHONIOENCODING=utf-8
run() {   # SET EXCL MODE [mep_reprocess args]
  SET=$1; EXCL=$2; MODE=$3; shift 3; export MEPSET=$SET
  PLOTS=""; [ "$MODE" = "full" ] || PLOTS="--no-plots"
  cd ${data_dir}/Revision_2026/Manuscript_Figures/code/mep_reprocess || exit 1
  python mep_reprocess.py "$R\extract" "$R\\$SET" --no-pulse "$R\no_pulse.csv" --intensity-overrides "$R\intensity_overrides.csv" --exclude-subjects "$EXCL" $PLOTS "$@" > "$RU/${SET}_console.txt" 2>&1; echo "$SET mep_reprocess $?"
  python "$RU/qc/make_table_mep_probability.py" "$SET" > "$RU/$SET/table_mep_probability_console.txt" 2>&1; echo "prob_table $?"
  cd "$C8" || exit 1
  A="$RU/$SET/analysis"; mkdir -p "$A"
  python mep_stats_A_vs_B_reprocessed.py > "$RU/$SET/analysis_mep_stats_console.txt" 2>&1; echo "stats $?"
  python tost_eo_ec_reprocessed.py > "$A/tost_console.txt" 2>&1; echo "tost $?"
  python build_single_trial_inputs_reprocessed.py > "$A/single_trial_inputs_console.txt" 2>&1; echo "st_inputs $?"
  if [ "$MODE" = "full" ]; then
    "$RS" Fig2_MEP_boxplot_reprocessed.R > "$A/fig2_console.txt" 2>&1; echo "fig2 $?"
    for f in single_trial_robustness_reprocessed.R fig6_single_trial_lmm_reprocessed.R fig6_single_trial_supp_reprocessed.R prestim_power_lmm_reprocessed.R prestim_power_lmm_v2_reprocessed.R mu_3e_reprocessed.R; do
      "$RS" $f > "$A/${f%.R}_console.txt" 2>&1; echo "$f $?"; done
    python mu_quintiles_reprocessed.py > "$A/mu_quintiles_console.txt" 2>&1; echo "quint $?"
  else
    for f in single_trial_robustness_reprocessed.R prestim_power_lmm_v2_reprocessed.R mu_3e_reprocessed.R; do
      "$RS" $f > "$A/${f%.R}_console.txt" 2>&1; echo "$f $?"; done
  fi
  WA="${data_dir}/Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/$SET/analysis"
  python fig6_input_reprocessed.py > /dev/null 2>&1; echo "fig6_input $?"
  "$RS" "${data_dir}/Revision_2026/Manuscript_Figures/code/Fig6_new_TEP_MEP.R" "$WA/Fig6_TEP_MEP_input_reprocessed.csv" "$WA/FigS16_TEP_MEP_reprocessed" "$WA/Fig6_correlations_reprocessed.csv" > "$A/figS16_console.txt" 2>&1; echo "S16 $?"
  "$RS" bootstrap_itpc_mep_17_reprocessed.R > "$A/bootstrap_console.txt" 2>&1; echo "boot $?"
  python corr_spec_ladder_reprocessed.py > "$A/ladder_console.txt" 2>&1; echo "ladder $?"
  python tep_mep_sliding_reprocessed.py > "$A/sliding_console.txt" 2>&1; echo "sliding $?"
  python fig6_pearson_logmep_reprocessed.py > "$A/pearson_console.txt" 2>&1; echo "pearson $?"
  (cd ${data_dir}/Revision_2026/Manuscript_Figures/code && python "$C8/eyestate_mep_rows_reprocessed.py" > "$A/eyestate_mep_rows_console.txt" 2>&1; echo "eyestate_rows $?")
  if [ "$MODE" = "full" ]; then
    python fig6_reprocessed.py > "$A/fig6_console.txt" 2>&1; echo "fig6 $?"
    python fig7_reprocessed.py > "$A/fig7_console.txt" 2>&1; echo "fig7 $?"
  fi
  echo "done $SET"
}
run main_qc2_008 014,025,040 full --qc-log "$V/qc_log_v2_main.csv"
run main_008 014,025,040 reduced
run main_qc2_absent_missing_008 014,025,040 reduced --qc-log "$V/qc_log_v2_absent_missing.csv"
run excl_artifact_qc2_008 014,025,040,006,018 reduced --qc-log "$V/qc_log_v2_main.csv"
cd ${data_dir}/Revision_2026/Manuscript_Figures
for s in main_qc2_008 main_qc2_absent_missing_008; do MEPSET=$s "$RS" "$W/late_window_mep_covariate_reprocessed.R" > "$RU/$s/analysis/late_window_console.txt" 2>&1; echo "late $s $?"; done
MEPSET=main_qc2_008 "$RS" "$C8/fig6_single_trial_combined_reprocessed.R" > "$RU/main_qc2_008/analysis/combined_console.txt" 2>&1; echo "combined $?"
cd "$RU/code"; export MEPSET=main_qc2_008
python "$C8/fig7_a3.py" > "$V/fig7_main_qc2_008_console.txt" 2>&1; echo "fig7_a3 $?"
python "$W/single_trial_per_participant_plot_reprocessed.py" > "$V/S18_main_qc2_008_console.txt" 2>&1; echo "S18 $?"
echo "F8 done"
