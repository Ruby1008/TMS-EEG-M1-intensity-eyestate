# qc/v2/make_code008.py copy (part 8): spTMS_008_Real_EC no longer excluded.
"""Part 2 step 7-4 (user 2026-10-22): single-trial inputs with the reprocessed MEP.
EEG trial <-> MEP trial: as the original programs (data/EEG_MEP_trial_map.csv, key B = trial2023 -> VETA trial_number). The MEP columns of key B are
replaced by the reprocessed values of <MEPSET>/mep_trials_reprocessed.csv: mep_int_B = intensity (after intensity_overrides.csv), mep_amp_B = p2p (mV),
mep_accept_B = accepted (1/0), new column mep_present_B = present (p2p >= 0.05 mV). Every accepted trial has p2p > 0, so no trial is dropped for MEP = 0.
Outputs (in <MEPSET>/analysis/): EEG_MEP_trial_map_reprocessed.csv, Single_trial_robustness_input_reprocessed.csv (data/Single_trial_robustness_input.csv
with the MEP columns replaced). Original files are read only."""
data_dir = r"path/to/data"   # release: set this path
import os
from pathlib import Path
import numpy as np, pandas as pd
MF = Path(data_dir + r"\Revision_2026\Manuscript_Figures"); SET = os.environ.get("MEPSET", "main")
RD = MF / "MEP_reprocess_2026-10-22" / SET; OUT = RD / "analysis"; OUT.mkdir(exist_ok=True)
R = pd.read_csv(RD / "mep_trials_reprocessed.csv", dtype={"subject": str})
R["sub"] = R.subject.astype(int)
key = R.set_index(["sub", "eye", "trial_number"])[["intensity", "p2p", "accepted", "present"]]
assert key.index.is_unique
M = pd.read_csv(MF / "MEP_reprocess_2026-10-22" / "qc" / "EEG_MEP_trial_map_C2fix.csv")   # part 4 A1: 022 EO key B corrected (qc/a1_022_mapping_check.py); previous version *_before_022fix.py
k = pd.MultiIndex.from_arrays([M["sub"], M.eye, M.trial2023])
M["mep_int_B"] = key.intensity.reindex(k).values; M["mep_amp_B"] = key.p2p.reindex(k).values
M["mep_accept_B"] = key.accepted.reindex(k).astype(float).values; M["mep_present_B"] = key.present.reindex(k).astype(float).values
assert not ((M.mep_accept_B == 1) & ~(M.mep_amp_B > 0)).any(), "accepted trial with p2p <= 0"
M.to_csv(OUT / "EEG_MEP_trial_map_reprocessed.csv", index=False)
X = pd.read_csv(MF / "data" / "Single_trial_robustness_input.csv").drop(columns=["mep_int_B", "mep_amp_B", "mep_accept_B"])
X = X.merge(M[["id", "epoch", "mep_int_B", "mep_amp_B", "mep_accept_B", "mep_present_B"]], on=["id", "epoch"], how="left")
X.to_csv(OUT / "Single_trial_robustness_input_reprocessed.csv", index=False)
use = M[M.mep_int_B.notna() & (M.mep_int_B == M.intensity) & (M.mep_accept_B == 1)]
print(SET, "| map rows", len(M), "| usable (intensity match, accepted, 008 EC excluded):", len(use), "| per intensity:", use.groupby("intensity").size().to_dict(),
      "| present share:", use.groupby("intensity").mep_present_B.mean().round(3).to_dict(), "| 022 EO intensity mismatches:", int(((M["sub"] == 22) & (M.eye == "EO") & (M.mep_int_B != M.intensity)).sum()))
