"""Section 3-1 follow-up (stop rule triggered, user 2026-10-12): where does the ERM pre-stimulus offset come from?
Variants, all on the Real main set (stage-D excluded removed), ROI FC3/FC1/C3/C1, per dataset x intensity:
  orig   ERM (IMF4-7 sum), Phase A shift definition (ref -300..-100, shift -20..-5, post 15..50)
  a      ERM + baseline correction -100..-10 ms (shift and post measured after that baseline)
  b      ERM + 20-ms sliding mean (21-sample moving average, i.e. the smoothing used for Fig 5) - the metrics
         in the main check were computed on the UNsmoothed ERM, so 'without smoothing' == orig; this shows the
         smoothed version instead
  bp     control: conventional _F_ERM (no low-pass) band-passed 4.2-36.4 Hz (4th-order Butterworth, zero phase)
         = the IMF4-7 band without EMD -> is the offset a property of band-limiting rather than of CEEMDAN?
  tep    conventional _F_TEP (reference)
  imf4..imf7   contribution of each IMF to the ERM shift (which mode carries the pre-stimulus lobe)
Tests per variant: ERM-vs-TEP |shift| (Wilcoxon), rmcorr(shift, post) over subjects, 12-cell shift==0 (t, Holm),
TOST +-0.5 / +-1 (Holm), plus the regression slope of shift on post.
Outputs: data/S1_variants_per_cell.csv, data/S1_variants_summary.csv, data/S1_variants_cells.csv, FigS_ERMvsTEP_variants.png
"""
from pathlib import Path
import numpy as np, pandas as pd, scipy.io as sio
from scipy import stats, signal
import pingouin as pg
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import os, sys as _s; _s.path.insert(0, str(Path(__file__).resolve().parent)); from tep_io import MAIN17
SUB17 = os.environ.get("SUBSET") == "main17"; SFX = "_main17" if SUB17 else ""   # 2026-10-16: 17-participant main version

ROOT = Path(__file__).resolve().parent.parent; HH = ROOT / "HHSA"; OUT = ROOT / "data"
Q = pd.read_csv(ROOT.parent / "Preproc" / "QC" / "QC_summary_long.csv", encoding="utf-8-sig")[["id", "intensity", "stageD"]]
W = pd.read_csv(OUT / "TEP_roi_waveforms.csv"); tcols = [c for c in W.columns if c.startswith("t_")]; tt = np.arange(-2000, -2000 + len(tcols))
LEV = [80, 100, 120]

def metrics(t, w, ref=(-300, -100)):
    w = w - w[(t >= ref[0]) & (t <= ref[1])].mean()
    return w[(t >= -20) & (t <= -5)].mean(), w[(t >= 15) & (t <= 50)].mean()

b, a = signal.butter(4, [4.2, 36.4], btype="band", fs=1000)
rows, waves = [], {}
for f in sorted(HH.glob("erm_spTMS_*_Real_*.mat")):
    R = sio.loadmat(f, squeeze_me=True); id_ = str(R["id"]); t = R["time"]; erm = R["erm_roi"]; imf = R["erm_roi_imf"]   # 3 x 4 x T
    for k, lv in enumerate(LEV):
        w = erm[k]; v = {}
        v["orig"] = metrics(t, w)
        v["a"] = metrics(t, w, ref=(-100, -10))
        ws = np.convolve(w, np.ones(21) / 21, mode="same"); v["b"] = metrics(t, ws)
        for j, im in enumerate([4, 5, 6, 7]): v[f"imf{im}"] = metrics(t, imf[k, j])
        g = W[(W.id == id_) & (W.intensity == lv)]
        if len(g):
            x_erm = g[g.file == "ERM"][tcols].to_numpy().ravel(); x_tep = g[g.file == "TEP"][tcols].to_numpy().ravel()
            v["bp"] = metrics(tt, signal.filtfilt(b, a, x_erm)); v["tep"] = metrics(tt, x_tep)
            if Q[(Q.id == id_) & (Q.intensity == lv)].stageD.iloc[0] != "exclude" and (not SUB17 or int(id_[6:9]) in MAIN17): waves.setdefault((id_[-2:], lv), []).append((w - w[(t >= -300) & (t <= -100)].mean(), signal.filtfilt(b, a, x_erm) - 0, x_tep))
        for name, (sh, po) in v.items():
            rows.append({"id": id_, "intensity": lv, "variant": name, "shift": sh, "post": po})
D = pd.DataFrame(rows); D["sub"] = D.id.str.extract(r"spTMS_(\d+)_").astype(int); D["eyes"] = D.id.str[-2:]
D = D.merge(Q, on=["id", "intensity"], how="left"); D = D[D.stageD != "exclude"]
if SUB17: D = D[D["sub"].isin(MAIN17)]
D.to_csv(OUT / f"S1_variants_per_cell{SFX}.csv", index=False)

def holm(p): return pg.multicomp(np.asarray(p, float), method="holm")[1]
def tost_p(x, bd):
    n = len(x); m = x.mean(); se = x.std(ddof=1) / np.sqrt(n)
    return max(stats.t.sf((m + bd) / se, n - 1), stats.t.cdf((m - bd) / se, n - 1))
tep = D[D.variant == "tep"].set_index(["id", "intensity"])["shift"]
summ, cells = [], []
for name in ["orig", "a", "b", "bp", "tep", "imf4", "imf5", "imf6", "imf7"]:
    d = D[D.variant == name].copy()
    r = pg.rm_corr(data=d, x="shift", y="post", subject="sub"); ci = r["CI95" if "CI95" in r else "CI95%"].iloc[0]
    sl = np.polyfit(d.post - d.groupby("sub").post.transform("mean"), d["shift"] - d.groupby("sub")["shift"].transform("mean"), 1)[0]
    cmp_ = d.set_index(["id", "intensity"])["shift"].abs() - tep.abs().reindex(d.set_index(["id", "intensity"]).index)
    pw = stats.wilcoxon(cmp_.dropna()).pvalue if name != "tep" else np.nan
    cell = []
    for e in ("EO", "EC"):
        for lv in LEV:
            x = d[(d.eyes == e) & (d.intensity == lv)]["shift"].to_numpy(); n = len(x); m = x.mean(); se = x.std(ddof=1) / np.sqrt(n)
            cell.append({"variant": name, "eyes": e, "intensity": lv, "n": n, "mean": m, "ci90_lo": m - stats.t.ppf(.95, n - 1) * se, "ci90_hi": m + stats.t.ppf(.95, n - 1) * se,
                         "p_t": stats.ttest_1samp(x, 0).pvalue, "p_tost_0.5": tost_p(x, .5), "p_tost_1": tost_p(x, 1)})
    C = pd.DataFrame(cell)
    for c in ("p_t", "p_tost_0.5", "p_tost_1"): C[c + "_holm"] = holm(C[c])
    cells.append(C)
    summ.append({"variant": name, "n_cells": len(d), "abs_shift_median_uV": d["shift"].abs().median(), "shift_sd_uV": d["shift"].std(),
                 "larger_than_TEP_pct": (cmp_ > 0).mean() * 100 if name != "tep" else np.nan, "p_vs_TEP_wilcoxon": pw,
                 "rmcorr_r": r["r"].iloc[0], "rmcorr_ci_lo": ci[0], "rmcorr_ci_hi": ci[1], "rmcorr_p": r["pval"].iloc[0],
                 "slope_shift_per_uV_post": sl, "cells_shift_ne0_holm": int((C.p_t_holm < .05).sum()), "cells_TOST05_holm": int((C["p_tost_0.5_holm"] < .05).sum()),
                 "cells_TOST1_holm": int((C.p_tost_1_holm < .05).sum())})
S = pd.DataFrame(summ); S.to_csv(OUT / f"S1_variants_summary{SFX}.csv", index=False); pd.concat(cells).to_csv(OUT / f"S1_variants_cells{SFX}.csv", index=False)
print(S.round(3).to_string())

# --- grand-average overlays: ERM orig / ERM+baseline / band-passed _F_ERM / _F_TEP, Real EO/EC x intensity
fig, axs = plt.subplots(2, 3, figsize=(14, 6.5), sharex=True)
for r_, e in enumerate(("EO", "EC")):
    for c_, lv in enumerate(LEV):
        ax = axs[r_, c_]; L = waves.get((e, lv), [])
        if not L: continue
        erm = np.mean([x[0] for x in L], 0); bp = np.mean([x[1] for x in L], 0); tp = np.mean([x[2] for x in L], 0)
        t = sio.loadmat(next(HH.glob("erm_spTMS_*_Real_*.mat")), squeeze_me=True)["time"]
        h1, = ax.plot(t, erm, color="#c55505", lw=2)
        h2, = ax.plot(t, erm - erm[(t >= -100) & (t <= -10)].mean(), color="#f9954e", lw=1.5, ls="--")
        h3, = ax.plot(tt, bp - bp[(tt >= -300) & (tt <= -100)].mean(), color="#2a78d6", lw=1.5)
        h4, = ax.plot(tt, tp - tp[(tt >= -300) & (tt <= -100)].mean(), color="k", lw=1.2)
        ax.legend([h1, h2, h3, h4], ["ERM IMF4-7 (ref −300..−100)", "(a) ERM, baseline −100..−10", "_F_ERM band-pass 4.2–36.4 Hz (no EMD)", "_F_TEP (30 Hz LP)"], fontsize=7, loc="upper right")
        ax.axvspan(-5, 10, color="k", alpha=.75); ax.axvspan(-20, -5, color="#ddd"); ax.axhline(0, color="#999", lw=.7)
        ax.set_xlim(-150, 300); ax.set_title(f"Real {e} {lv}% RMT")
axs[1, 0].set_xlabel("Time (ms)"); axs[0, 0].set_ylabel("ROI (µV)"); axs[1, 0].set_ylabel("ROI (µV)")
fig.suptitle("Grand averages (ROI); grey band = pre-stimulus window −20..−5 ms", fontsize=11); fig.tight_layout()
fig.savefig(ROOT / f"FigS_ERMvsTEP_variants{SFX}.png", dpi=600)
