"""Supplementary TEP measures (master plan C; decisions 2026-10-16), for ERM, _F_TEP (TEP) and the band-pass control (BP).
Per recording x intensity (all available real recordings; data/ALLCH_avg.npz and HHSA/erm_*.mat):
  fixed      P30/N45/P60/N100 mean in the grand-average windows (data/TEP_component_windows.csv)
  indiv      P30/N45/P60: each participant's own peak within the literature range (P30 20-40, N45 38-55, P60 56-82 ms) in the average of that
             participant's available conditions, mean +-10 ms around it per condition [definition chosen here: one latency per participant]
  gmfp       mean global mean field power 15-60 ms (61 scalp channels; SD across channels at each sample)
  area       ROI area 15-60 ms (sum of |ROI| per ms, uV*ms)
rmANOVA (2 x 3, GG; pingouin) on the 17 main participants for indiv/gmfp/area (fixed: tep_stats.py). LMM on all available participants: R
(tep_supp_lmm.R). Outputs: data/TEP_supp_measures_long.csv, data/TEP_supp_indiv_latencies.csv, data/TEP_supp_anova.csv
"""
import numpy as np, pandas as pd, pingouin as pg
import tep_io as io

Z = np.load(io.DATA / "ALLCH_avg.npz"); labs = list(Z["labels"]); t = np.round(Z["time"]).astype(int)
scalp = [i for i, l in enumerate(labs) if l not in ("M1", "M2")]; ri = [labs.index(c) for c in io.ROI]
W = pd.read_csv(io.DATA / "TEP_component_windows.csv")
RANGES = {"P30": (20, 40, 1), "N45": (38, 55, -1), "P60": (56, 82, 1)}
def data(kind, s, e):
    """(time, ROI chan-avg [time x 3], scalp [chan x time x 3])"""
    if kind == "ERM":
        R = io.load_erm(s, e); te = np.round(R["time"]).astype(int)
        return te, R["erm_roi"].T, np.transpose(R["erm_allch"], (0, 2, 1))     # erm_allch: chan x 3 x time -> chan x time x 3
    A = Z[f"spTMS_{s:03d}_Real_{e}|{kind}"].astype(float); return t, A[ri].mean(0), A[scalp]
rows, lat = [], []
for kind in ("ERM", "TEP", "BP"):
    win = {r.component: (r.win_start_ms, r.win_end_ms) for r in W[W.data == kind].itertuples()}
    for s in sorted(set(io.AVAIL["EO"]) | set(io.AVAIL["EC"])):
        eyes = [e for e in io.EYES if s in io.AVAIL[e]]
        D = {e: data(kind, s, e) for e in eyes}; tt = D[eyes[0]][0]
        pav = np.mean([D[e][1][:, k] for e in eyes for k in range(3)], 0)       # participant average over available conditions
        ilat = {}
        for c, (a, b, sg) in RANGES.items():
            w = np.where((tt >= a) & (tt <= b))[0]; ilat[c] = int(tt[w[np.argmax(sg * pav[w])]])
            lat.append(dict(data=kind, sub=s, component=c, latency_ms=ilat[c], edge=ilat[c] in (a, b)))
        for e in eyes:
            tt, roi, sc = D[e]
            for k, lv in enumerate(io.LEV):
                base = dict(data=kind, sub=s, eye=e, intensity=lv, main17=s in io.MAIN17)
                for c, (a0, a1) in win.items(): rows.append(dict(base, measure=f"fixed_{c}", value=roi[(tt >= a0) & (tt <= a1), k].mean()))
                for c, L in ilat.items(): rows.append(dict(base, measure=f"indiv_{c}", value=roi[(tt >= L - 10) & (tt <= L + 10), k].mean()))
                m = (tt >= 15) & (tt <= 60)
                rows.append(dict(base, measure="gmfp_15_60", value=sc[:, m, k].std(0).mean()))
                rows.append(dict(base, measure="area_15_60", value=np.abs(roi[m, k]).sum()))
L = pd.DataFrame(rows); L.to_csv(io.DATA / "TEP_supp_measures_long.csv", index=False)
LT = pd.DataFrame(lat); LT.to_csv(io.DATA / "TEP_supp_indiv_latencies.csv", index=False)
print("individual peaks on a range edge:", LT.groupby(["data", "component"]).edge.sum().to_dict())
A = []
for (kind, meas), g in L[L.main17].groupby(["data", "measure"]):
    if meas.startswith("fixed"): continue
    a = pg.rm_anova(data=g, dv="value", within=["eye", "intensity"], subject="sub", correction=True, effsize="np2").set_index("Source")
    r = dict(data=kind, measure=meas, n=g["sub"].nunique())
    for src, k in (("eye", "eye"), ("intensity", "intensity"), ("eye * intensity", "interaction")):
        x = a.loc[src]; r.update({f"F_{k}": x["F"], f"df1_{k}": x["ddof1"], f"df2_{k}": x["ddof2"], f"p_{k}": x["p_GG_corr"] if pd.notna(x.get("p_GG_corr", np.nan)) else x["p_unc"],
                                  f"np2_{k}": x["np2"], f"eps_{k}": x["eps"] if pd.notna(x.get("eps", np.nan)) else 1.0})
    A.append(r)
A = pd.DataFrame(A); A.to_csv(io.DATA / "TEP_supp_anova.csv", index=False)
pd.set_option("display.width", 220); print(A[["data", "measure", "n", "F_intensity", "eps_intensity", "p_intensity", "np2_intensity", "p_eye", "p_interaction"]].round(4).to_string())
