"""Manipulation check (user 2026-10-18, section 3.1; pre-specified analysis added for the reviewers): pre-stimulus alpha power, EO vs EC.
Data: 17 main participants; real stimulation (Preproc\\Final\\<id>_F_ERM.set, 0.1-80 Hz, average reference) and sham (same participants where available).
Window -1000..-10 ms (991 samples) of every accepted trial; regions: occipital O1/O2/POz, sensorimotor ROI FC3/FC1/C3/C1.
 (1) main: Welch PSD per trial and channel (Hann, 500-sample segments, 50% overlap, 2-Hz bins), mean 8-12 Hz (bins 8, 10, 12 Hz) -> mean over channels of the
     region and over trials -> log10 (µV²/Hz);
 (2) IMF6 version: mean squared first-layer CEEMDAN IMF6 (HHSA\\timf_spm_<id>.dat, same window) -> mean over channels and trials -> log10 (µV²).
Intensities pooled (main; intensities are randomly interleaved within a block) and per intensity (table).
Statistics: paired t-test EC vs EO (difference EC - EO), Cohen's dz = mean / SD of the paired difference; Holm over the 2 pooled tests (occipital, sensorimotor) within each stimulation type (main, user 2026-10-19) and over the 4 pooled tests (kept).
Outputs: data/Alpha_prestim_values.csv, data/Alpha_prestim_tests.csv, FigS_alpha_prestim.png/.tiff (600 dpi)."""
import numpy as np, pandas as pd, scipy.io as sio
from scipy import stats, signal
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import tep_io as io

REG = {"occipital": ["O1", "O2", "POz"], "sensorimotor": io.ROI}
def timf_imf6(id_, chans):
    S = sio.loadmat(io.HHSA / f"timf_spm_{id_}.mat", squeeze_me=True, struct_as_record=False)["D"]
    dim = tuple(int(v) for v in S.data.dim); assert int(S.data.dtype) == 16 and dim[1] == 11
    labs = [c.label for c in np.atleast_1d(S.channels)]; t = (S.timeOnset + np.arange(int(S.Nsamples)) / float(S.Fsample)) * 1000
    w = np.where((np.round(t) >= -1000) & (np.round(t) <= -10))[0]; ci = [labs.index(c) for c in chans]
    nch, nimf, nt, ntr = dim; blk = nch * nimf; x = np.empty((len(ci), len(w), ntr))
    with open(io.HHSA / f"timf_spm_{id_}.dat", "rb") as fh:                                   # one contiguous read per trial (Fortran order: chan, imf, time, trial)
        for k in range(ntr):
            fh.seek(4 * (blk * nt * k + blk * w[0])); a = np.fromfile(fh, dtype="<f4", count=blk * len(w)).reshape((nch, nimf, len(w)), order="F")
            x[:, :, k] = a[ci, 5, :]                                                                  # IMF6
    
    inten = np.array([float(tr.label) for tr in np.atleast_1d(S.trials)]); return x, inten
rows, spec = [], {}
for stim in ("Real", "Sham"):
    for s in io.MAIN17:
        for e in io.EYES:
            id_ = f"spTMS_{s:03d}_{stim}_{e}"
            if not (io.FINAL / f"{id_}_F_ERM.set").exists(): continue
            x, labs, t, inten, _ = io.load_set(io.FINAL, id_ + "_F_ERM"); t = np.round(t).astype(int)
            w = (t >= -1000) & (t <= -10); assert w.sum() == 991
            f, P = signal.welch(x[:, w, :].astype(float), fs=1000, window="hann", nperseg=500, noverlap=250, axis=1)   # chan x freq x trial
            fa = (f >= 8) & (f <= 12)
            try: xi, inten_i = timf_imf6(id_, sum(REG.values(), []))
            except FileNotFoundError: xi = None
            print(id_, flush=True)
            for reg, ch in REG.items():
                ci = [labs.index(c) for c in ch]; pw = P[ci][:, fa].mean((0, 1))                              # per trial
                spec[(stim, s, e, reg)] = P[ci].mean((0, 2))
                if xi is not None:
                    k0 = sum(REG.values(), []).index(ch[0]); pi = (xi[k0:k0 + len(ch)] ** 2).mean((0, 1)); assert np.allclose(np.sort(inten_i), np.sort(inten))
                for lv in ("all", 80, 100, 120):
                    m = np.ones(len(inten), bool) if lv == "all" else inten == lv
                    rows.append(dict(stim=stim, sub=s, eye=e, region=reg, intensity=str(lv), method="Welch", n_trials=int(m.sum()), power=pw[m].mean(), log10_power=np.log10(pw[m].mean())))
                    if xi is not None:
                        mi = np.ones(len(inten_i), bool) if lv == "all" else inten_i == lv
                        rows.append(dict(stim=stim, sub=s, eye=e, region=reg, intensity=str(lv), method="IMF6", n_trials=int(mi.sum()), power=pi[mi].mean(), log10_power=np.log10(pi[mi].mean())))
V = pd.DataFrame(rows); V.to_csv(io.DATA / "Alpha_prestim_values.csv", index=False)
tests = []
for (stim, reg, lv, meth), g in V.groupby(["stim", "region", "intensity", "method"]):
    w_ = g.pivot(index="sub", columns="eye", values="log10_power").dropna(); d = w_["EC"] - w_["EO"]; r = stats.ttest_rel(w_["EC"], w_["EO"])
    tests.append(dict(stim=stim, region=reg, intensity=str(lv), method=meth, n=len(d), mean_EO=w_["EO"].mean(), mean_EC=w_["EC"].mean(), mean_diff_EC_minus_EO=d.mean(),
                      t=r.statistic, df=len(d) - 1, p=r.pvalue, dz=d.mean() / d.std(ddof=1), n_EC_greater=int((d > 0).sum())))
T = pd.DataFrame(tests); T["p_holm_pooled"] = np.nan
for meth in ("Welch", "IMF6"):
    k = (T.intensity == "all") & (T.method == meth)
    if k.any():
        p = T.loc[k, "p"].values; o = np.argsort(p); h = np.minimum(1, np.maximum.accumulate(p[o] * (len(p) - np.arange(len(p))))); hh = np.empty_like(h); hh[o] = h; T.loc[k, "p_holm_pooled"] = hh
# 2026-10-19 (user decision 2-c): main correction = Holm over the 2 pooled tests (occipital, sensorimotor) within each stimulation type; the 4-test version is kept
T["p_holm_within_stim"] = np.nan
for meth in ("Welch", "IMF6"):
    for stim in ("Real", "Sham"):
        k = (T.intensity == "all") & (T.method == meth) & (T.stim == stim)
        if k.any():
            p = T.loc[k, "p"].values; o = np.argsort(p); h = np.minimum(1, np.maximum.accumulate(p[o] * (len(p) - np.arange(len(p))))); hh = np.empty_like(h); hh[o] = h; T.loc[k, "p_holm_within_stim"] = hh
T.to_csv(io.DATA / "Alpha_prestim_tests.csv", index=False)
print(T[T.intensity == "all"].round(4).to_string(index=False))
# ---- figure: spectra (real, sham) and paired values (Welch, pooled)
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 7}); C = {"EO": "#E08A2C", "EC": "#3B8C8C"}
fig, ax = plt.subplots(2, 4, figsize=(7.0, 3.6))
for i, stim in enumerate(("Real", "Sham")):
    for j, reg in enumerate(REG):
        a = ax[i, j]
        for e in io.EYES:
            S_ = np.array([v for (st, s, ee, rg), v in spec.items() if st == stim and ee == e and rg == reg]); fm = (f >= 2) & (f <= 40)
            a.semilogy(f[fm], S_.mean(0)[fm], color=C[e], lw=1.2, label=e)
        a.axvspan(8, 12, color="0.92", zorder=0); a.set_xlabel("Frequency (Hz)"); a.set_title(f"{stim}, {reg}"); a.legend(frameon=False, fontsize=6)
        if j == 0: a.set_ylabel("PSD (µV²/Hz), −1000 to −10 ms")
        b = ax[i, j + 2]; g = V[(V.stim == stim) & (V.region == reg) & (V.intensity == "all") & (V.method == "Welch")].pivot(index="sub", columns="eye", values="log10_power").dropna()
        for _, r in g.iterrows(): b.plot([0, 1], [r.EO, r.EC], color="0.7", lw=.6, marker="o", ms=2.5, mfc="0.4", mec="none")
        b.plot([0, 1], [g.EO.mean(), g.EC.mean()], color="k", lw=1.5, marker="o", ms=4)
        tt = T[(T.stim == stim) & (T.region == reg) & (T.intensity == "all") & (T.method == "Welch")].iloc[0]
        b.set_xticks([0, 1]); b.set_xticklabels(["EO", "EC"]); b.set_xlim(-.4, 1.4); b.set_title(f"{stim}, {reg}", fontsize=7); b.text(.5, .97, "***" if tt.p_holm_within_stim < .001 else "**" if tt.p_holm_within_stim < .01 else "*" if tt.p_holm_within_stim < .05 else "", transform=b.transAxes, ha="center", va="top", fontsize=11)
        if j == 0: b.set_ylabel("log10 alpha power (8–12 Hz)")
fig.tight_layout()
for ext in ("png", "tiff"): fig.savefig(io.ROOT / f"FigS_alpha_prestim.{ext}", dpi=600, **({"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tiff" else {}))
print("saved")
