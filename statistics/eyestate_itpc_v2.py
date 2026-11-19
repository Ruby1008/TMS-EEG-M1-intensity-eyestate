"""Round 15 merge (2026-10-22): ITPC rows of the eye-state analyses (eyestate_null_round8.py sections 3c (i), (iv), 3f and eyestate_null_round8.R
3c (ii), (v)) recomputed with ITPC rule v2 (baseline -500..-100 ms; windows = longest significant intensity segment of rule v2, data/ITPC_v2_windows.csv).
Same functions and models as the originals (copied); only the ITPC values differ. Outputs: data/EyeState_allavailable_long_ITPC_v2.csv,
data/EyeState_null_robustness_ITPC_v2.csv ((i), (iv); (ii), (v) appended by eyestate_itpc_v2.R), data/EyeState_individual_differences_v2.csv."""
import subprocess
import numpy as np, pandas as pd, pingouin as pg, scipy.io as sio
from scipy import stats
import tep_io as io

D = io.DATA; S17 = io.MAIN17; ROWS = []
def add(**k): ROWS.append(k)
def rm(df, dv, within=("eye", "intensity")):
    a = pg.rm_anova(data=df, dv=dv, within=list(within), subject="sub", correction=True, effsize="np2").set_index("Source"); return a
def gg(a, src): r = a.loc[src]; return r["F"], r["ddof1"], r["ddof2"], (r["p_GG_corr"] if "p_GG_corr" in r and pd.notna(r.get("p_GG_corr")) else r["p_unc"]), r["np2"]
M = sio.loadmat(D / "ITPC_roi_real.mat", squeeze_me=True); ids = [str(i) for i in M["ids"]]; t = np.round(M["t"].astype(float)); idx = {i: k for k, i in enumerate(ids)}
X = M["X"] - M["X"][:, :, (t >= -500) & (t <= -100), :].mean(2, keepdims=True)
WN = pd.read_csv(D / "ITPC_v2_windows.csv").set_index("IMF"); IW = {imf: (WN.loc[imf, "longest_start_ms"], WN.loc[imf, "longest_end_ms"]) for imf in (6, 7)}
long = []
for e in io.EYES:
    for s in sorted(set(io.AVAIL["EO"]) | set(io.AVAIL["EC"])):
        if s not in io.AVAIL[e]: continue
        k = idx[f"spTMS_{s:03d}_Real_{e}"]
        for li, lv in enumerate(io.LEV):
            for imf in (6, 7): long.append(dict(sub=s, eye=e, intensity=lv, measure=f"ITPC IMF{imf}", value=X[k, imf - 1][(t >= IW[imf][0]) & (t <= IW[imf][1]), li].mean()))
L = pd.DataFrame(long); L["main17"] = L["sub"].isin(S17); L.to_csv(D / "EyeState_allavailable_long_ITPC_v2.csv", index=False)
# check: main-sample values equal data/Fig8_ITPC_N100_input_v2.csv
F8 = pd.read_csv(D / "Fig8_ITPC_N100_input_v2.csv").set_index("sub")
for r in L[L.main17].itertuples(): assert np.isclose(r.value, F8.loc[r.sub, f"ITPC_{r.measure.split()[1]}_{r.eye}_{r.intensity}"])
_q = io.Q[(io.Q.status == "processed") & (io.Q.stim == "Real") & (io.Q.stageD == "sensitivity")]; SENS = [s for s in S17 if s not in set(_q["sub"])]
for label, subs_eeg in (("(i) main analysis", S17), ("(iv) sensitivity, no 20-50 uV residual", SENS)):
    for meas in ["ITPC IMF6", "ITPC IMF7"]:
        g = L[(L.measure == meas) & L["sub"].isin(subs_eeg)]; full = g.groupby("sub").size(); g = g[g["sub"].isin(full[full == 6].index)]
        a = rm(g, "value")
        for src in ("eye", "eye * intensity"):
            F, d1, d2, p, np2 = gg(a, src); add(section="3c", analysis=label, measure=meas, effect=src, n=g["sub"].nunique(), F=F, df1=d1, df2=d2, p=p, np2=np2)
pd.DataFrame(ROWS).to_csv(D / "EyeState_null_robustness_ITPC_v2.csv", index=False)
subprocess.run(["Rscript", str(io.ROOT / "code" / "eyestate_itpc_v2.R")], check=True)
# 3f individual differences: IMF6 ITPC row recomputed; log MEP and ERM N100 rows from data/EyeState_individual_differences.csv; Holm over the three
LAPd = pd.read_csv(D / "EyeState_laplacian_values.csv")
lap = LAPd[(LAPd.stim == "Real") & (LAPd.version == "after Laplacian (CSD)") & (LAPd.region == "C3")].pivot(index="sub", columns="eye", values="log10_power"); dmu = lap.EC - lap.EO
g = L[(L.measure == "ITPC IMF6") & L["sub"].isin(S17)].pivot_table(index=["sub", "intensity"], columns="eye", values="value"); y = (g.EC - g.EO).groupby("sub").mean()
j = dmu.index.intersection(y.dropna().index); r = stats.spearmanr(dmu.loc[j], y.loc[j])
ID0 = pd.read_csv(D / "EyeState_individual_differences.csv"); ID = ID0[ID0.y != "EC - EO ITPC IMF6"][["x", "y", "n", "rho", "p"]].copy()
ID = pd.concat([ID, pd.DataFrame([dict(x="EC - EO Laplacian mu (C3)", y="EC - EO ITPC IMF6 (rule v2)", n=len(j), rho=r.statistic, p=r.pvalue)])], ignore_index=True)
ID["p_holm"] = pg.multicomp(ID.p.values, method="holm")[1]; ID.to_csv(D / "EyeState_individual_differences_v2.csv", index=False)
pd.set_option("display.width", 200); print(pd.read_csv(D / "EyeState_null_robustness_ITPC_v2.csv").round(4).to_string()); print(ID.round(3).to_string())
