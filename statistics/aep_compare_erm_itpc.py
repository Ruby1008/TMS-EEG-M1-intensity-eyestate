"""AEP-ICA robustness (master plan J, S3-S5): ERM and ITPC before vs after AEP removal, 17 participants.
Before = HHSA\\erm_*.mat (main analysis); after = HHSA_AEPr\\erm_*.mat (aep_remove.m -> run_layer1_aepr.m, same CEEMDAN parameters/seed).
ERM: sliding-window 2 x 3 rmANOVA (21 samples, 10-ms steps, 0-400 ms, GG, BH-FDR per effect) and N100 (120-140 ms) ANOVA, before and after.
ITPC: Fig 7 analysis (itpc_fig78_G.py; after = ITPC_VARIANT=AEPr) -> overlap of the FDR masks per effect and IMF.
P30/N45/P60: after the component windows are decided.
Outputs: data/AEP_ERM_slidingwindow_before_after.csv, data/AEP_ERM_N100_anova.csv, data/AEP_ITPC_mask_overlap.csv,
         FigS3_ERM_AEP_before_after.png/.tiff (600 dpi)
"""
import numpy as np, pandas as pd, pingouin as pg
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import tep_io as io

S17, LEV, EYES = io.MAIN17, io.LEV, io.EYES
FOLD = {"before": io.HHSA, "after": io.ROOT / "HHSA_AEPr"}
ERM = {v: {(s, e): io.load_erm(s, e, f) for s in S17 for e in EYES} for v, f in FOLD.items()}
t = np.round(ERM["before"][(S17[0], "EO")]["time"]).astype(int); starts = np.arange(0, 400, 10)
def bh(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); q = np.empty(m); q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]; return np.minimum(q, 1)
def anova(df):
    a = pg.rm_anova(data=df, dv="y", within=["eye", "intensity"], subject="sub", correction=True, effsize="np2").set_index("Source"); out = {}
    for src, k in (("eye", "eye"), ("intensity", "intensity"), ("eye * intensity", "interaction")):
        r = a.loc[src]; out.update({f"F_{k}": r["F"], f"p_{k}": r["p_GG_corr"] if pd.notna(r.get("p_GG_corr", np.nan)) else r["p_unc"], f"np2_{k}": r["np2"],
                                    f"eps_{k}": r["eps"] if pd.notna(r.get("eps", np.nan)) else 1.0})
    return out
SW, N1 = [], []
for v in FOLD:
    rows = []
    for a0 in starts:
        i0 = np.where(t >= a0)[0][0]
        df = pd.DataFrame([dict(sub=s, eye=e, intensity=lv, y=ERM[v][(s, e)]["erm_roi"][k][i0:i0 + 21].mean()) for s in S17 for e in EYES for k, lv in enumerate(LEV)])
        rows.append(dict(version=v, win_start_ms=a0, **anova(df)))
    R = pd.DataFrame(rows)
    for k in ("eye", "intensity", "interaction"): R[f"q_{k}"] = bh(R[f"p_{k}"])
    SW.append(R)
    w = (t >= 120) & (t <= 140)
    df = pd.DataFrame([dict(sub=s, eye=e, intensity=lv, y=ERM[v][(s, e)]["erm_roi"][k][w].mean()) for s in S17 for e in EYES for k, lv in enumerate(LEV)])
    N1.append(dict(version=v, component="N100", window_ms="120-140", n=len(S17), **anova(df)))
SW = pd.concat(SW); SW.to_csv(io.DATA / "AEP_ERM_slidingwindow_before_after.csv", index=False)
N1 = pd.DataFrame(N1); N1.to_csv(io.DATA / "AEP_ERM_N100_anova.csv", index=False)
def runs(m):
    out, st = [], None
    for a0, x in zip(starts, m):
        if x and st is None: st = a0
        if not x and st is not None: out.append(f"{st}-{prev + 20}"); st = None
        prev = a0
    if st is not None: out.append(f"{st}-{prev + 20}")
    return ", ".join(out) or "none"
for v in FOLD:
    r = SW[SW.version == v]; print(v, "ERM intensity FDR:", runs(r.q_intensity < .05), "| eye FDR:", runs(r.q_eye < .05), "| interaction FDR:", runs(r.q_interaction < .05))
pd.set_option("display.width", 220); print(N1.round(4).to_string())
# ITPC masks
B = pd.read_csv(io.DATA / "Fig7_ITPC_anova_maps.csv"); A = pd.read_csv(io.DATA / "Fig7_ITPC_anova_maps_AEPr.csv")
M = B.merge(A, on=["effect", "IMF", "time_ms"], suffixes=("_before", "_after"))
ov = []
for (ef, imf), g in M.groupby(["effect", "IMF"]):
    b, a = g.qBH_before < .05, g.qBH_after < .05
    ov.append(dict(effect=ef, IMF=imf, sig_before=int(b.sum()), sig_after=int(a.sum()), both=int((a & b).sum()), before_only=int((b & ~a).sum()), after_only=int((a & ~b).sum()),
                   jaccard=(a & b).sum() / max((a | b).sum(), 1), r_F=np.corrcoef(g.F_before, g.F_after)[0, 1]))
OV = pd.DataFrame(ov); OV.to_csv(io.DATA / "AEP_ITPC_mask_overlap.csv", index=False); print(OV.round(3).to_string())
# Fig S3
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 9})
fig, axs = plt.subplots(2, 3, figsize=(11, 5.8), sharex=True, sharey=True)
for r_, e in enumerate(EYES):
    for c_, lv in enumerate(LEV):
        ax = axs[r_, c_]; k = LEV.index(lv)
        for v, sty in (("before", dict(color="k", lw=1.3)), ("after", dict(color="#c0392b", lw=1.3, ls="--"))):
            ax.plot(t, np.mean([ERM[v][(s, e)]["erm_roi"][k] for s in S17], 0), label=f"{'Before' if v == 'before' else 'After'} AEP removal", **sty)
        for v, y, col in (("before", 5.5, "k"), ("after", 5.0, "#c0392b")):
            for a0 in starts[(SW[SW.version == v].q_intensity < .05).values]: ax.plot([a0, a0 + 20], [y, y], color=col, lw=3, solid_capstyle="butt")
        ax.axvspan(120, 140, color="0.92", lw=0); ax.axhline(0, color="0.6", lw=.5); ax.set_xlim(0, 400); ax.set_title(f"Real {e}, {lv}% RMT")
        if c_ == 0: ax.set_ylabel("ERM, ROI (µV)")
        if r_ == 1: ax.set_xlabel("Time (ms)")
axs[0, 0].legend(frameon=False, fontsize=7, loc="lower right")
fig.suptitle("ERM (IMF4–7) before and after AEP-ICA removal", fontsize=9)
fig.tight_layout(); fig.savefig(io.ROOT / "FigS3_ERM_AEP_before_after.png", dpi=600); fig.savefig(io.ROOT / "FigS3_ERM_AEP_before_after.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
