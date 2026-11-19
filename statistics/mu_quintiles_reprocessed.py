data_dir = r"path/to/data"   # release: set this path
# qc/v2/make_code008.py copy (part 8): spTMS_008_Real_EC no longer excluded.
# copy made by MEP_reprocess_2026-10-22/code/make_reprocessed_copies.py (2026-10-22); input = reprocessed MEP data (MEPSET = main | excl_artifact); outputs in MEP_reprocess_2026-10-22/<MEPSET>/analysis/.
import sys; sys.path.insert(0, data_dir + r'\Revision_2026\Manuscript_Figures\code')
"""Round 9 (user 2026-10-21, section 3, Figure 7 D–F). Real stimulation, 17 participants.
D: for every scalp channel, single-trial log10 8–12 Hz power (−1000..−10 ms; data/Prestim_power_single_trial_allch.npz) centred within participant × intensity
   (power_within), ROI ERM P30 (27–47 ms; data/ERM_single_trial.csv) ~ power_within + intensity + (1 | participant), all epochs, maximum likelihood (statsmodels MixedLM);
   t of power_within; BH-FDR over 61 channels -> data/SingleTrial_mu_channel_P30.csv
E: ROI mu power (data/Prestim_power_single_trial.csv, log_mu) in within-participant quintiles per intensity; P30 centred within participant × intensity; mean per
   participant and quintile, then mean ± SE over participants -> data/SingleTrial_mu_quintiles.csv (outcome = P30)
F: as E with log MEP (data/Single_trial_robustness_input.csv: accepted trials with a discernible MEP and matching intensity, 008 EC excluded; 100% and 120% RMT,
   as the single-trial models of Figure 7 A–C), log MEP centred within participant × intensity -> same file (outcome = log MEP)"""
import numpy as np, pandas as pd, warnings
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import tep_io as io
import os as _os
from pathlib import Path as _P
_SET = _os.environ.get("MEPSET", "main"); _RD = _P(data_dir + r"\Revision_2026\Manuscript_Figures\MEP_reprocess_2026-10-22") / _SET
_OUT = _RD / "analysis"; _OUT.mkdir(exist_ok=True)
warnings.filterwarnings("ignore")
D = io.DATA
Z = np.load(D / "Prestim_power_single_trial_allch.npz"); CH = list(Z["channels"]); ids = list(Z["ids"])
ST = pd.read_csv(D / "ERM_single_trial.csv")
rows = []
for i in ids:
    P = Z[f"P_{i}"]; I = Z[f"I_{i}"]; st = ST[ST.id == i].sort_values("epoch"); assert len(st) == len(I) and np.array_equal(st.intensity.values, I)
    d = pd.DataFrame(P, columns=CH); d["P30"] = st.P30.values; d["intensity"] = I; d["sub"] = int(i.split("_")[1]); d["id"] = i; rows.append(d)
A = pd.concat(rows, ignore_index=True)
O = pd.read_csv(D / "SingleTrial_mu_channel_P30.csv")   # D: P30 side only, unchanged (not recomputed)
out = []
for c in []:
    x = A[["sub", "intensity", "P30", c]].rename(columns={c: "pw"}).copy()
    x["pw_w"] = x.pw - x.groupby(["sub", "intensity"]).pw.transform("mean")
    m = smf.mixedlm("P30 ~ pw_w + C(intensity)", x, groups=x["sub"]).fit(reml=False)
    out.append(dict(channel=c, b=m.params["pw_w"], se=m.bse["pw_w"], t=m.tvalues["pw_w"], p=m.pvalues["pw_w"], n_trials=len(x), sd_pw_w=x.pw_w.std()))
pass
print("D: channels q < .05:", int((O.q_fdr61 < .05).sum()), "; ROI t:", O[O.channel.isin(io.ROI)].t.round(2).tolist())
# E / F quintiles
PW = pd.read_csv(D / "Prestim_power_single_trial.csv"); X = pd.read_csv(_OUT / "Single_trial_robustness_input_reprocessed.csv")
E = PW.merge(ST[["id", "epoch", "P30"]], on=["id", "epoch"]); E["sub"] = E.id.str.split("_").str[1].astype(int)
def quint(df, y):
    df = df.copy(); df["q"] = df.groupby(["sub", "intensity"]).log_mu.transform(lambda v: pd.qcut(v.rank(method="first"), 5, labels=False) + 1)
    df["yc"] = df[y] - df.groupby(["sub", "intensity"])[y].transform("mean")
    per = df.groupby(["intensity", "sub", "q"]).yc.mean().reset_index()
    return per.groupby(["intensity", "q"]).yc.agg(mean="mean", se=lambda v: v.std(ddof=1) / np.sqrt(len(v)), n="count").reset_index()
qE = quint(E, "P30").assign(outcome="P30 (µV, centred)")
F = X[X.mep_int_B.notna() & (X.mep_int_B == X.intensity) & (X.mep_accept_B == 1) & (X.mep_amp_B > 0) & X.intensity.isin([100, 120])].copy()
F["logMEP"] = np.log(F.mep_amp_B); qF = quint(F, "logMEP").assign(outcome="log MEP (centred)")
Q = pd.concat([qE, qF]); Q.to_csv(_OUT / "SingleTrial_mu_quintiles_reprocessed.csv", index=False); print(Q.round(3).to_string(index=False))
