"""R2-1 reliable-response analysis (master plan I), real vs sham, the 17 main participants (sham: same participants where available; 024 sham EO missing).
Data: Preproc\\Final\\<id>_F_TEP.set (0.1-30 Hz), 61 scalp channels (M1/M2 excluded); data version and all details below accepted by the user 2026-10-17.
Per recording x intensity:
  (1) GMFP test: observed = mean GMFP (SD across channels of the trial average) over 15-60 ms (46 samples). Null: 2000 bootstrap draws, each trial contributes a
      random 46-ms segment starting in -1000..-56 ms (ending <= -10 ms), the segments are averaged over the same trials and the mean GMFP is computed;
      p = (1 + #null >= observed) / 2001.
  (2) single-trial RMS (over channels and samples) 15-60 ms vs -100..-55 ms, one-sided Wilcoxon signed-rank (post > pre).
  Benjamini-Hochberg FDR over all recording x intensity cells (real and sham together) separately for each test; reliable response = both q < .05.
  log SNR = log(observed GMFP / mean of the null distribution).
Group tests: responder counts per stimulation x eye x intensity; Cochran's Q across intensities (participants with all three cells) per stimulation x eye;
McNemar (exact, discordant pairs) real vs sham per eye x intensity; log SNR real vs sham per eye x intensity, paired Wilcoxon, Holm across six.
Outputs: data/R21_cells.csv, data/R21_group_tests.csv, FigS_R21_GA_real_sham.png/.tiff, FigS_R21_responders.png/.tiff, FigS_R21_logSNR.png/.tiff (600 dpi)
"""
import numpy as np, pandas as pd, pingouin as pg
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import tep_io as io

import os
# Robustness variant (user 2026-10-18, section 4.1): R21_SHAM_VARIANT=noOverride -> sham recordings read from Revision_2026\Preproc_Sham_noOverride\Final
# when present there (the 21 sham recordings with manual ICA decisions; the others are identical to the main pipeline); outputs get the suffix _shamNoOverride.
# Real cells are processed first with the same seed, so they are identical to the main analysis. Unset = main analysis (unchanged).
VAR = os.environ.get("R21_SHAM_VARIANT", ""); SFX = "_shamNoOverride" if VAR == "noOverride" else ""; NOOVR = io.REV / "Preproc_Sham_noOverride" / "Final"
RNG = np.random.default_rng(20261017); B = 2000; L = 46
rows, GA = [], {}
for stim in ("Real", "Sham"):
    for s in io.MAIN17:
        for e in io.EYES:
            id_ = f"spTMS_{s:03d}_{stim}_{e}"
            if not (io.FINAL / f"{id_}_F_TEP.set").exists(): continue
            fold = NOOVR if (SFX and stim == "Sham" and (NOOVR / f"{id_}_F_TEP.set").exists()) else io.FINAL
            x, labs, t, inten, _ = io.load_set(fold, id_ + "_F_TEP"); t = np.round(t).astype(int)
            sc = [i for i, l in enumerate(labs) if l not in ("M1", "M2")]; ri = [labs.index(c) for c in io.ROI]; x = x[sc].astype(float)
            post = np.where((t >= 15) & (t <= 60))[0]; assert len(post) == L
            starts = np.where((t >= -1000) & (t <= -56))[0]; pre_rms = np.where((t >= -100) & (t <= -55))[0]
            for lv in io.LEV:
                X = x[:, :, inten == lv]; n = X.shape[2]
                obs = X[:, post, :].mean(2).std(0).mean()
                Xt = np.transpose(X, (2, 1, 0))                               # trial x time x chan
                null = np.empty(B)
                for b0 in range(0, B, 100):
                    st = RNG.choice(starts, size=(100, n))                     # draws x trials
                    idx = st[:, :, None] + np.arange(L)                        # draws x trials x L
                    seg = Xt[np.arange(n)[None, :, None], idx]                 # draws x trials x L x chan
                    null[b0:b0 + 100] = seg.mean(1).std(2).mean(1)
                p_boot = (1 + (null >= obs).sum()) / (B + 1)
                rms_post = np.sqrt((X[:, post, :] ** 2).mean((0, 1))); rms_pre = np.sqrt((X[:, pre_rms, :] ** 2).mean((0, 1)))
                p_w = stats.wilcoxon(rms_post, rms_pre, alternative="greater").pvalue
                rows.append(dict(id=id_, sub=s, stim=stim, eye=e, intensity=lv, n_trials=n, gmfp_post=obs, gmfp_null_mean=null.mean(), log_snr=np.log(obs / null.mean()),
                                 p_boot=p_boot, p_rms=p_w, source=fold.parent.name))
                GA.setdefault((stim, e, lv), []).append(X[ri].mean(0).mean(1))
            print(id_)
C = pd.DataFrame(rows)
def bh(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); q = np.empty(m); q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]; return np.minimum(q, 1)
C["q_boot"] = bh(C.p_boot); C["q_rms"] = bh(C.p_rms); C["responder"] = (C.q_boot < .05) & (C.q_rms < .05)
C.to_csv(io.DATA / f"R21_cells{SFX}.csv", index=False)
G = []
for (stim, e, lv), g in C.groupby(["stim", "eye", "intensity"]): G.append(dict(test="responders", stim=stim, eye=e, intensity=lv, n=len(g), k=int(g.responder.sum()), value=g.responder.mean()))
for (stim, e), g in C.groupby(["stim", "eye"]):
    W = g.pivot(index="sub", columns="intensity", values="responder").dropna().astype(int)
    k = W.sum(1); cj = W.sum(0); N = W.values.sum(); kk = 3
    Q = (kk - 1) * (kk * (cj ** 2).sum() - N ** 2) / (kk * N - (k ** 2).sum()) if (kk * N - (k ** 2).sum()) > 0 else np.nan
    G.append(dict(test="Cochran Q (intensity)", stim=stim, eye=e, n=len(W), value=Q, p=stats.chi2.sf(Q, kk - 1) if np.isfinite(Q) else np.nan))
mc = []
for e in io.EYES:
    for lv in io.LEV:
        a = C[(C.stim == "Real") & (C.eye == e) & (C.intensity == lv)].set_index("sub"); b = C[(C.stim == "Sham") & (C.eye == e) & (C.intensity == lv)].set_index("sub")
        sub = a.index.intersection(b.index); ra, rb = a.loc[sub, "responder"], b.loc[sub, "responder"]
        n10, n01 = int((ra & ~rb).sum()), int((~ra & rb).sum())
        p_mc = stats.binomtest(n10, n10 + n01, .5).pvalue if n10 + n01 > 0 else 1.0
        w = stats.wilcoxon(a.loc[sub, "log_snr"], b.loc[sub, "log_snr"])
        mc.append(dict(eye=e, intensity=lv, n_pairs=len(sub), real_only=n10, sham_only=n01, p_mcnemar=p_mc, logsnr_real=a.loc[sub, "log_snr"].mean(), logsnr_sham=b.loc[sub, "log_snr"].mean(), W=w.statistic, p_logsnr=w.pvalue))
MC = pd.DataFrame(mc); MC["p_logsnr_holm"] = pg.multicomp(MC.p_logsnr.values, method="holm")[1]
for r in MC.itertuples():
    G.append(dict(test="McNemar real vs sham", eye=r.eye, intensity=r.intensity, n=r.n_pairs, value=f"{r.real_only}/{r.sham_only}", p=r.p_mcnemar))
    G.append(dict(test="log SNR real vs sham (Wilcoxon)", eye=r.eye, intensity=r.intensity, n=r.n_pairs, value=f"{r.logsnr_real:.3f} vs {r.logsnr_sham:.3f}", p=r.p_logsnr, p_holm=r.p_logsnr_holm))
GT = pd.DataFrame(G); GT.to_csv(io.DATA / f"R21_group_tests{SFX}.csv", index=False)
if SFX: raise SystemExit("variant: tables written, figures skipped")
pd.set_option("display.width", 200); print(GT.to_string())
# ---- figures
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 9}); t = np.arange(-2000, 2000)
fig, axs = plt.subplots(2, 3, figsize=(11, 5.6), sharex=True, sharey=True)
for r_, e in enumerate(io.EYES):
    for c_, lv in enumerate(io.LEV):
        ax = axs[r_, c_]
        for stim, col in (("Real", "#c0392b"), ("Sham", "0.3")):
            L_ = GA.get((stim, e, lv), [])
            if L_: ax.plot(t, np.mean(L_, 0), color=col, lw=1.3, label=stim)
        ax.axvspan(-5, 10, color="k", lw=0); ax.axvspan(15, 60, color="0.92", lw=0, zorder=0); ax.axhline(0, color="0.6", lw=.5); ax.set_xlim(-100, 400); ax.set_title(f"{e}, {lv}% RMT"); ax.legend(frameon=False, fontsize=7)
        if c_ == 0: ax.set_ylabel("ROI amplitude (µV)")
        if r_ == 1: ax.set_xlabel("Time (ms)")
fig.suptitle("Conventional TEP (0.1–30 Hz), ROI FC3/FC1/C3/C1: real vs sham stimulation; grey: 15–60 ms", fontsize=10); fig.tight_layout()
fig.savefig(io.ROOT / "FigS_R21_GA_real_sham.png", dpi=600); fig.savefig(io.ROOT / "FigS_R21_GA_real_sham.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"}); plt.close(fig)
R = GT[GT.test == "responders"]
fig, axs = plt.subplots(1, 2, figsize=(8, 3.4), sharey=True)
for ax, e in zip(axs, io.EYES):
    for k_, (stim, col) in enumerate((("Real", "#c0392b"), ("Sham", "0.5"))):
        g = R[(R.stim == stim) & (R.eye == e)].sort_values("intensity"); xs = np.arange(3) + (k_ - .5) * .38
        ax.bar(xs, g.value, width=.36, color=col, label=stim)
        for xx, (kk_, nn) in zip(xs, zip(g.k, g.n)): ax.text(xx, 0.02 + (kk_ / nn), f"{int(kk_)}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(range(3)); ax.set_xticklabels(["80%", "100%", "120%"]); ax.set_title(e); ax.set_ylim(0, 1.1)
axs[0].set_ylabel("Proportion of reliable responses"); axs[0].legend(frameon=False); fig.tight_layout()
fig.savefig(io.ROOT / "FigS_R21_responders.png", dpi=600); fig.savefig(io.ROOT / "FigS_R21_responders.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"}); plt.close(fig)
fig, axs = plt.subplots(1, 2, figsize=(8, 3.4), sharey=True)
for ax, e in zip(axs, io.EYES):
    for stim, col in (("Real", "#c0392b"), ("Sham", "0.4")):
        g = C[(C.stim == stim) & (C.eye == e)].groupby("intensity").log_snr; m, se = g.mean(), g.sem()
        ax.errorbar(m.index, m.values, yerr=se.values, color=col, marker="o", capsize=3, label=stim)
    ax.axhline(0, color="0.6", lw=.5); ax.set_xticks(io.LEV); ax.set_xlabel("Intensity (% RMT)"); ax.set_title(e)
axs[0].set_ylabel("log SNR (mean ± SEM)"); axs[0].legend(frameon=False); fig.tight_layout()
fig.savefig(io.ROOT / "FigS_R21_logSNR.png", dpi=600); fig.savefig(io.ROOT / "FigS_R21_logSNR.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"}); plt.close(fig)
