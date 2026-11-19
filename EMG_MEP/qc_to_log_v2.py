"""Part 5 E3 (user 2026-10-23): qc/qc_to_log.py extended for QC v2 (decision unsure). The v1 file is not changed.
Usage: python qc_to_log_v2.py <exported.csv> <out_dir> [--map ../qc_id_map.csv] [--prefix qc_log_v2]
Checks: every qc_id of qc_id_map.csv exactly once, decision in {accept, absent, reject, unsure}, a reason for every absent / reject / unsure.
Outputs (subject, eye, trial_number, decision, reason; accept -> empty decision and reason; input format of mep_reprocess.py --qc-log):
  <prefix>_main.csv              unsure -> accept (main analysis)
  <prefix>_unsure_absent.csv     unsure -> absent (reason "unsure-as-absent: <reason>")
  <prefix>_absent_missing.csv    absent and unsure -> reject (reason "absent-as-missing: ..." / "unsure-as-missing: ...")
Prints the counts of the four decisions (blinded: no decoding of subject, eye or intensity here)."""
import argparse, sys
from pathlib import Path
import pandas as pd
ap = argparse.ArgumentParser(); ap.add_argument("exported"); ap.add_argument("out_dir")
ap.add_argument("--map", default=str(Path(__file__).parent.parent / "qc_id_map.csv")); ap.add_argument("--prefix", default="qc_log_v2")
a = ap.parse_args()
E = pd.read_csv(a.exported, dtype=str, encoding="utf-8-sig", keep_default_na=False); M = pd.read_csv(a.map, dtype={"subject": str})
if list(E.columns[:3]) != ["qc_id", "decision", "reason"]: sys.exit(f"columns must be qc_id, decision, reason (found {list(E.columns)})")
E["qc_id"] = E.qc_id.str.strip(); E["decision"] = E.decision.str.strip().str.lower(); E["reason"] = E.reason.str.strip()
err = []
dup = E.qc_id[E.qc_id.duplicated()].unique()
if len(dup): err.append(f"duplicated qc_id: {', '.join(dup[:20])}")
unknown = sorted(set(E.qc_id) - set(M.qc_id)); missing = sorted(set(M.qc_id) - set(E.qc_id))
if unknown: err.append(f"{len(unknown)} unknown qc_id: {', '.join(unknown[:20])}")
if missing: err.append(f"{len(missing)} qc_id missing: {', '.join(missing[:20])}")
bad = E[~E.decision.isin(["accept", "absent", "reject", "unsure"])]
if len(bad): err.append(f"{len(bad)} invalid decisions: {', '.join(bad.qc_id[:20])}")
nor = E[E.decision.isin(["absent", "reject", "unsure"]) & E.reason.eq("")]
if len(nor): err.append(f"{len(nor)} absent/reject/unsure without a reason: {', '.join(nor.qc_id[:20])}")
if err: sys.exit("QC CSV not converted:\n  " + "\n  ".join(err))
L = M.merge(E[["qc_id", "decision", "reason"]], on="qc_id", how="left"); od = Path(a.out_dir); od.mkdir(parents=True, exist_ok=True)
cols = ["subject", "eye", "trial_number", "decision", "reason"]
def save(df, name): df = df.copy(); df.loc[df.decision.eq("accept"), ["decision", "reason"]] = ""; df[cols].to_csv(od / f"{a.prefix}_{name}.csv", index=False)
m = L.copy(); u = m.decision.eq("unsure"); m.loc[u, "decision"] = "accept"; save(m, "main")
x = L.copy(); u = x.decision.eq("unsure"); x.loc[u, "decision"] = "absent"; x.loc[u, "reason"] = "unsure-as-absent: " + x.loc[u, "reason"]; save(x, "unsure_absent")
y = L.copy()
for d in ("absent", "unsure"):
    k = y.decision.eq(d); y.loc[k, "reason"] = f"{d}-as-missing: " + y.loc[k, "reason"]; y.loc[k, "decision"] = "reject"
save(y, "absent_missing")
c = E.decision.value_counts()
print(f"converted {len(E)} trials: accept {c.get('accept', 0)}, absent {c.get('absent', 0)}, reject {c.get('reject', 0)}, unsure {c.get('unsure', 0)} -> {od}\\{a.prefix}_main / _unsure_absent / _absent_missing.csv")
