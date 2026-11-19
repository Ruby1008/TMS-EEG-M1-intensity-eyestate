#!/usr/bin/env bash
data_dir="path/to/data"   # release: set this path
# Part 4 (user 2026-10-22): the part-2 pipeline (mep_reprocess.py step 5/6 + step 7 analysis copies) in the order used on 2026-10-22, for one data set.
# Usage: bash run_pipeline.sh <SET> <exclude-subjects> [full|reduced] [extra mep_reprocess.py arguments ...]
#   full    = all step-7 programs and figures (as main); reduced = the programs run for excl_artifact on 2026-10-22 (+ figures off)
# Example: bash run_pipeline.sh main 014,025,040 full
#          bash run_pipeline.sh main_qc 014,025,040 full --qc-log "<data_dir>\...\qc\qc_log_user_2026-10-29.csv"
set -u
SET=$1; EXCL=$2; MODE=${3:-full}; shift 3 || true
R="${data_dir}\Revision_2026\Manuscript_Figures\MEP_reprocess_2026-10-22"
RU="${data_dir}/Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22"
RS="Rscript"
export PYTHONIOENCODING=utf-8 MEPSET=$SET
PLOTS=""; [ "$MODE" = "full" ] || PLOTS="--no-plots"
cd ${data_dir}/Revision_2026/Manuscript_Figures/code/mep_reprocess || exit 1
python mep_reprocess.py "$R\extract" "$R\\$SET" --no-pulse "$R\no_pulse.csv" --intensity-overrides "$R\intensity_overrides.csv" --exclude-subjects "$EXCL" $PLOTS "$@" > "$RU/${SET}_console.txt" 2>&1; echo "mep_reprocess $?"
python "$RU/qc/make_table_mep_probability.py" "$SET" > "$RU/$SET/table_mep_probability_console.txt" 2>&1; echo "prob_table $?"
cd "$RU/code" || exit 1
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
"$RS" ../../code/Fig6_new_TEP_MEP.R "$WA/Fig6_TEP_MEP_input_reprocessed.csv" "$WA/FigS16_TEP_MEP_reprocessed" "$WA/Fig6_correlations_reprocessed.csv" > "$A/figS16_console.txt" 2>&1; echo "S16 $?"
"$RS" bootstrap_itpc_mep_17_reprocessed.R > "$A/bootstrap_console.txt" 2>&1; echo "boot $?"
python corr_spec_ladder_reprocessed.py > "$A/ladder_console.txt" 2>&1; echo "ladder $?"
python tep_mep_sliding_reprocessed.py > "$A/sliding_console.txt" 2>&1; echo "sliding $?"
python fig6_pearson_logmep_reprocessed.py > "$A/pearson_console.txt" 2>&1; echo "pearson $?"
(cd ${data_dir}/Revision_2026/Manuscript_Figures/code && python "$RU/code/eyestate_mep_rows_reprocessed.py" > "$A/eyestate_mep_rows_console.txt" 2>&1; echo "eyestate_rows $?")
if [ "$MODE" = "full" ]; then
  python fig6_reprocessed.py > "$A/fig6_console.txt" 2>&1; echo "fig6 $?"
  python fig7_reprocessed.py > "$A/fig7_console.txt" 2>&1; echo "fig7 $?"
fi
echo "done $SET"
