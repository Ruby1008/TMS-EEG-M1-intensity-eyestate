#!/usr/bin/env bash
data_dir="path/to/data"   # release: set this path
# Part 5 F1: QC v2 versions (run_pipeline.sh of part 4; outputs in MEP_reprocess_2026-10-22/<SET>/), late-window and combined models, figures of main_qc2.
Q=${data_dir}/Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc
V="${data_dir}/Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc/v2"
W="${data_dir}/Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22/qc"
RS="Rscript"
cd "$Q"
bash run_pipeline.sh main_qc2 014,025,040 full --qc-log "$V/qc_log_v2_main.csv"
bash run_pipeline.sh main_qc2_unsure_absent 014,025,040 reduced --qc-log "$V/qc_log_v2_unsure_absent.csv"
bash run_pipeline.sh main_qc2_absent_missing 014,025,040 reduced --qc-log "$V/qc_log_v2_absent_missing.csv"
bash run_pipeline.sh excl_artifact_qc2 014,025,040,006,018 reduced --qc-log "$V/qc_log_v2_main.csv"
bash run_pipeline.sh main_qc2_018seg7 014,025,040 reduced --qc-log "$V/qc_log_v2_main_018seg7.csv"
bash run_pipeline.sh main_qc2_018seg58 014,025,040 reduced --qc-log "$V/qc_log_v2_main_018seg58.csv"
cd ${data_dir}/Revision_2026/Manuscript_Figures
for s in main_qc2 main_qc2_unsure_absent main_qc2_absent_missing; do MEPSET=$s "$RS" "$W/late_window_mep_covariate_reprocessed.R" > "$Q/../$s/analysis/late_window_console.txt" 2>&1; echo "late $s $?"; done
for s in main_qc2 main_qc2_unsure_absent; do MEPSET=$s "$RS" "$W/fig6_single_trial_combined_reprocessed.R" > "$Q/../$s/analysis/combined_console.txt" 2>&1; echo "combined $s $?"; done
cd "$Q/../code"; export MEPSET=main_qc2 PYTHONIOENCODING=utf-8
"$RS" "$W/figures_a3/Fig2_MEP_boxplot_a3.R" > "$Q/v2/fig2_main_qc2_console.txt" 2>&1; echo "fig2_a3 $?"
python "$W/figures_a3/fig6_a3.py" > "$Q/v2/fig6_main_qc2_console.txt" 2>&1; echo "fig6_a3 $?"
python "$W/figures_a3/fig7_a3.py" > "$Q/v2/fig7_main_qc2_console.txt" 2>&1; echo "fig7_a3 $?"
python "$W/single_trial_per_participant_plot_reprocessed.py" > "$Q/v2/S18_main_qc2_console.txt" 2>&1; echo "S18 $?"
cd ${data_dir}/Revision_2026/Manuscript_Figures; MEPSET=main_qc2 "$RS" "$W/v2/Fig2_p6_qc2.R" > "$Q/v2/fig2_p6_main_qc2_console.txt" 2>&1; echo "fig2_p6 $?"
echo "F1 done"
