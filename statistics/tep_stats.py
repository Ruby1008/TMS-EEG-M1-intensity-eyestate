"""TEP statistics, master plan C (2026-10-16), 17 participants.
1. Sliding-window 2 (eye) x 3 (intensity) rmANOVA as spTMS_ERM.m: mean of 21 samples starting at 0, 10, ..., 390 ms (ERM: the
   erm_roi_win values of erm_extract.m; _F_TEP: same windows from data/TEP_roi_waveforms.csv); Greenhouse-Geisser p (pingouin);
   BH-FDR across the 40 windows per effect (main) and uncorrected (Supplementary). ERM = main (Fig 5B red line); _F_TEP and BP (_F_ERM band-passed 4.2-36.4 Hz, data/ALLCH_avg.npz) = S-E.
2. Component rmANOVA (supplementary to the sliding window): fixed windows from data/TEP_component_windows.csv when present
   (P30/N45/P60: grand-average peak +-10 ms, pending the user's decision); N100 = 120-140 ms (fixed); post hoc paired t on
   intensity (eye states pooled), Holm across the three contrasts.
Outputs: data/TEP_slidingwindow_anova_<ERM|TEP>.csv, data/TEP_component_values_<ERM|TEP>.csv, data/TEP_component_anova.csv,
         data/TEP_component_posthoc.csv
"""
import numpy as np, pandas as pd, pingouin as pg
from scipy import stats
import tep_io as io

LEV, EYES, S17 = io.LEV, io.EYES, io.MAIN17
import os; SET = os.environ.get('TEP_SET', '')   # 2026-10-16: 'sens' = sensitivity set (17 minus participants with any 20-50 uV condition)
if SET == 'sens':
    _q = io.Q[(io.Q.status == 'processed') & (io.Q.stim == 'Real') & (io.Q.stageD == 'sensitivity')]
    S17 = [s for s in S17 if s not in set(_q['sub'])]
SFX = '_sens' if SET == 'sens' else ''
starts = np.arange(0, 400, 10)
# ---- per participant x condition waveforms (ROI)
ERM = {(s, e): io.load_erm(s, e) for s in S17 for e in EYES}
t_erm = np.round(ERM[(S17[0], "EO")]["time"]).astype(int)
W = pd.read_csv(io.DATA / "TEP_roi_waveforms.csv"); W = W[W.file == "TEP"]; tc = [c for c in W.columns if c.startswith("t_")]; t_tep = np.arange(-2000, -2000 + len(tc))
TEP = {(int(r.id[6:9]), r.id[-2:], int(r.intensity)): np.asarray(r[tc], float) for _, r in W.iterrows() if r.id.split("_")[2] == "Real"}
Z = np.load(io.DATA / "ALLCH_avg.npz"); _l = list(Z["labels"]); _ri = [_l.index(c) for c in io.ROI]; t_bp = np.round(Z["time"]).astype(int)
BP = {(s, e): Z[f"spTMS_{s:03d}_Real_{e}|BP"][_ri].mean(0) for s in S17 for e in EYES}   # S-E control (2026-10-16): _F_ERM 4.2-36.4 Hz
def wave(kind, s, e, lv):
    if kind == "ERM": return (t_erm, ERM[(s, e)]["erm_roi"][LEV.index(lv)])
    if kind == "BP": return (t_bp, BP[(s, e)][:, LEV.index(lv)])
    return (t_tep, TEP[(s, e, lv)])
def winmean(t, y, a): i0 = np.where(t >= a)[0][0]; return y[i0:i0 + 21].mean()

def bh(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); q = np.empty(m)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]; return np.minimum(q, 1)
def anova(df, dv):
    a = pg.rm_anova(data=df, dv=dv, within=["eye", "intensity"], subject="sub", correction=True, effsize="np2").set_index("Source")
    out = {}
    for src, key in (("eye", "eye"), ("intensity", "intensity"), ("eye * intensity", "interaction")):
        r = a.loc[src]; pgg = r["p_GG_corr"] if "p_GG_corr" in r and pd.notna(r.get("p_GG_corr")) else r["p_unc"]
        out.update({f"F_{key}": r["F"], f"df1_{key}": r["ddof1"], f"df2_{key}": r["ddof2"], f"p_{key}": pgg, f"np2_{key}": r["np2"],
                    f"eps_{key}": r["eps"] if "eps" in r and pd.notna(r.get("eps")) else 1.0})
    return out

for kind in ("ERM", "TEP", "BP"):
    rows = []
    for a0 in starts:
        df = pd.DataFrame([dict(sub=s, eye=e, intensity=lv, y=winmean(*wave(kind, s, e, lv), a0)) for s in S17 for e in EYES for lv in LEV])
        rows.append(dict(win_start_ms=a0, win_end_ms=a0 + 20, **anova(df, "y")))
    R = pd.DataFrame(rows)
    for k in ("eye", "intensity", "interaction"): R[f"q_{k}"] = bh(R[f"p_{k}"])
    R.to_csv(io.DATA / f"TEP_slidingwindow_anova_{kind}{SFX}.csv", index=False)
    def runs(mask):
        out, st = [], None
        for a0, m in zip(starts, mask):
            if m and st is None: st = a0
            if not m and st is not None: out.append(f"{st}-{prev + 20}"); st = None
            prev = a0
        if st is not None: out.append(f"{st}-{prev + 20}")
        return ", ".join(out) or "none"
    print(kind, "| intensity FDR:", runs(R.q_intensity < .05), "| uncorrected:", runs(R.p_intensity < .05),
          "| eye FDR:", runs(R.q_eye < .05), "| interaction FDR:", runs(R.q_interaction < .05))

# ---- components
CW = io.DATA / "TEP_component_windows.csv"
comp = {"ERM": {"N100": (120, 140)}, "TEP": {"N100": (120, 140)}, "BP": {"N100": (120, 140)}}
if CW.exists():
    for r in pd.read_csv(CW).itertuples(): comp[r.data][r.component] = (r.win_start_ms, r.win_end_ms)
A, PH, VAL = [], [], []
for kind, cc in comp.items():
    for name, (a0, a1) in cc.items():
        df = pd.DataFrame([dict(sub=s, eye=e, intensity=lv, y=(lambda t, y: y[(t >= a0) & (t <= a1)].mean())(*wave(kind, s, e, lv))) for s in S17 for e in EYES for lv in LEV])
        VAL.append(df.assign(data=kind, component=name, window=f"{a0}-{a1}"))
        A.append(dict(data=kind, component=name, window_ms=f"{a0}-{a1}", n=len(S17), **anova(df, "y")))
        P = df.groupby(["sub", "intensity"]).y.mean().unstack()
        tmp = []
        for x, y in ((120, 100), (120, 80), (100, 80)):
            r = stats.ttest_rel(P[x], P[y]); tmp.append(dict(data=kind, component=name, contrast=f"{x} vs {y}", mean_a=P[x].mean(), mean_b=P[y].mean(), t=r.statistic, df=len(P) - 1, p=r.pvalue))
        ph = pd.DataFrame(tmp); ph["p_holm"] = pg.multicomp(ph.p.values, method="holm")[1]; PH.append(ph)
pd.concat(VAL).to_csv(io.DATA / f"TEP_component_values{SFX}.csv", index=False)
A = pd.DataFrame(A); A.to_csv(io.DATA / f"TEP_component_anova{SFX}.csv", index=False)
PH = pd.concat(PH); PH.to_csv(io.DATA / f"TEP_component_posthoc{SFX}.csv", index=False)
pd.set_option("display.width", 250)
print(A[["data", "component", "window_ms", "F_intensity", "df1_intensity", "df2_intensity", "eps_intensity", "p_intensity", "np2_intensity", "F_eye", "p_eye", "F_interaction", "p_interaction"]].round(4).to_string())
print(PH.round(4).to_string())
