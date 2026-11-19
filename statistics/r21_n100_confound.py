"""N100 confound (master plan I): real vs sham ERM (CEEMDAN IMF4-7, first layer, same parameters) at the ROI and at Cz, and a sliding-window 2 (eye) x 3 (intensity)
rmANOVA on the real - sham difference wave (21-sample windows, 10-ms steps, 0-400 ms, Greenhouse-Geisser, BH-FDR across windows per effect).
Participants: the 17 main participants with sham data in both eye states (n = 16; 024 sham EO missing). ERM files: HHSA\\erm_spTMS_<sub>_<Real|Sham>_<eye>.mat.
Also the N100 (120-140 ms) real vs sham per condition (paired t, Holm across six).
Outputs: data/R21_diffwave_slidingwindow.csv, data/R21_N100_real_vs_sham.csv, FigS_R21_ERM_real_sham_ROI.png/.tiff, FigS_R21_ERM_real_sham_Cz.png/.tiff (600 dpi)
"""
import numpy as np, pandas as pd, pingouin as pg, scipy.io as sio
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import tep_io as io

def erm(s, stim, e): return sio.loadmat(io.HHSA / f"erm_spTMS_{s:03d}_{stim}_{e}.mat", squeeze_me=True)
SUBS = [s for s in io.MAIN17 if all((io.HHSA / f"erm_spTMS_{s:03d}_Sham_{e}.mat").exists() for e in io.EYES)]
R0 = erm(SUBS[0], "Real", "EO"); t = np.round(R0["time"]).astype(int); labs = [str(x) for x in R0["chanlabels"]]; cz = labs.index("Cz")
W = {}
for s in SUBS:
    for stim in ("Real", "Sham"):
        for e in io.EYES:
            R = erm(s, stim, e); W[(s, stim, e, "ROI")] = R["erm_roi"]; W[(s, stim, e, "Cz")] = R["erm_allch"][cz]       # 3 x time
starts = np.arange(0, 400, 10)
def bh(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); q = np.empty(m); q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]; return np.minimum(q, 1)
rows = []
for a0 in starts:
    i0 = np.where(t >= a0)[0][0]
    df = pd.DataFrame([dict(sub=s, eye=e, intensity=lv, y=(W[(s, "Real", e, "ROI")][k] - W[(s, "Sham", e, "ROI")][k])[i0:i0 + 21].mean()) for s in SUBS for e in io.EYES for k, lv in enumerate(io.LEV)])
    a = pg.rm_anova(data=df, dv="y", within=["eye", "intensity"], subject="sub", correction=True, effsize="np2").set_index("Source"); r = dict(win_start_ms=a0)
    for src, k in (("eye", "eye"), ("intensity", "intensity"), ("eye * intensity", "interaction")):
        x = a.loc[src]; r.update({f"F_{k}": x["F"], f"p_{k}": x["p_GG_corr"] if pd.notna(x.get("p_GG_corr", np.nan)) else x["p_unc"], f"np2_{k}": x["np2"]})
    rows.append(r)
R = pd.DataFrame(rows)
for k in ("eye", "intensity", "interaction"): R[f"q_{k}"] = bh(R[f"p_{k}"])
R.to_csv(io.DATA / "R21_diffwave_slidingwindow.csv", index=False)
w = (t >= 120) & (t <= 140); N = []
for e in io.EYES:
    for k, lv in enumerate(io.LEV):
        a = np.array([W[(s, "Real", e, "ROI")][k][w].mean() for s in SUBS]); b = np.array([W[(s, "Sham", e, "ROI")][k][w].mean() for s in SUBS]); r = stats.ttest_rel(a, b)
        N.append(dict(eye=e, intensity=lv, n=len(SUBS), real=a.mean(), sham=b.mean(), t=r.statistic, p=r.pvalue))
N = pd.DataFrame(N); N["p_holm"] = pg.multicomp(N.p.values, method="holm")[1]; N.to_csv(io.DATA / "R21_N100_real_vs_sham.csv", index=False)
def runs(m):
    out, st, prev = [], None, None
    for a0, x in zip(starts, m):
        if x and st is None: st = a0
        if not x and st is not None: out.append(f"{st}-{prev + 20}"); st = None
        prev = a0
    if st is not None: out.append(f"{st}-{prev + 20}")
    return ", ".join(out) or "none"
print("n =", len(SUBS), "| diff-wave intensity FDR:", runs(R.q_intensity < .05), "| eye FDR:", runs(R.q_eye < .05), "| interaction FDR:", runs(R.q_interaction < .05))
pd.set_option("display.width", 200); print(N.round(4).to_string(index=False))
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 9})
for site in ("ROI", "Cz"):
    fig, axs = plt.subplots(2, 3, figsize=(11, 5.6), sharex=True, sharey=True)
    for r_, e in enumerate(io.EYES):
        for c_, lv in enumerate(io.LEV):
            ax = axs[r_, c_]; k = io.LEV.index(lv)
            for stim, col in (("Real", "#c0392b"), ("Sham", "0.3")): ax.plot(t, np.mean([W[(s, stim, e, site)][k] for s in SUBS], 0), color=col, lw=1.3, label=stim)
            ax.plot(t, np.mean([W[(s, "Real", e, site)][k] - W[(s, "Sham", e, site)][k] for s in SUBS], 0), color="#2a78d6", lw=1, ls="--", label="Real − Sham")
            if site == "ROI":
                for a0 in starts[(R.q_intensity < .05).values]: ax.plot([a0, a0 + 20], [ax.get_ylim()[1] if False else 4.5] * 2, color="#2a78d6", lw=3)
            ax.axvspan(120, 140, color="0.92", lw=0, zorder=0); ax.axhline(0, color="0.6", lw=.5); ax.set_xlim(0, 400); ax.set_title(f"{e}, {lv}% RMT")
            if c_ == 0: ax.set_ylabel(f"ERM, {site} (µV)")
            if r_ == 1: ax.set_xlabel("Time (ms)")
    axs[0, 0].legend(frameon=False, fontsize=7)
    fig.suptitle(f"ERM (IMF4–7) real vs sham stimulation, {site}", fontsize=9)
    fig.tight_layout(); fig.savefig(io.ROOT / f"FigS_R21_ERM_real_sham_{site}.png", dpi=600); fig.savefig(io.ROOT / f"FigS_R21_ERM_real_sham_{site}.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"}); plt.close(fig)
