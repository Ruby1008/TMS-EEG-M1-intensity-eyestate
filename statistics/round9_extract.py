"""Round 9 (user 2026-10-21, section 3): data for the new panels of Figures 6 and 7 and for Figure S1. Read-only on Preproc\\Final.
(1) Per-channel pre-stimulus alpha power (−1000..−10 ms, Welch 500-sample Hann, 50% overlap, 8–12 Hz, mean PSD over trials, log10), average reference and after the
    spherical-spline surface Laplacian (same CSD matrix as eyestate_null_round8.py: m = 4, 50 terms, λ = 1e-5, 63 channels), real and sham, 17 main participants;
    paired t EC − EO per channel (61 scalp channels), BH-FDR over 61 channels per stimulation × version.
    -> data/EyeState_channel_alpha.csv, data/EyeState_channel_alpha_t.csv
(2) Single-trial log10 mu-band (8–12 Hz) power per scalp channel (average reference), real stimulation, 17 participants, all epochs of _F_ERM (epoch = 1-based index,
    as data/Prestim_power_single_trial.csv) -> data/Prestim_power_single_trial_allch.npz
(3) Grand-average conventional TEP (_F_TEP, 0.1–30 Hz) per channel, eye state and intensity, 17 participants (participant means, then grand mean)
    -> data/TEP_GA_allch.npz (Figure S1)"""
data_dir = r"path/to/data"   # release: set this path
import numpy as np, pandas as pd
from scipy import signal, stats
from scipy.special import eval_legendre
import tep_io as io

D = io.DATA; S17 = io.MAIN17
CED = pd.read_csv(data_dir + r"\waveguard64_woRef.ced", sep="\t", skipinitialspace=True); CED.columns = [c.strip() for c in CED.columns]; CED["labels"] = CED["labels"].str.strip()
def csd_matrix(labs, m=4, nterms=50, lam=1e-5):           # identical to eyestate_null_round8.py
    P = np.array([[float(CED[CED.labels == l].X.iloc[0]), float(CED[CED.labels == l].Y.iloc[0]), float(CED[CED.labels == l].Z.iloc[0])] for l in labs]); P /= np.linalg.norm(P, axis=1, keepdims=True)
    cosang = np.clip(P @ P.T, -1, 1); n = np.arange(1, nterms + 1)
    Ln = np.stack([eval_legendre(k, cosang) for k in n])
    G = np.tensordot((2 * n + 1) / (n * (n + 1)) ** m, Ln, 1) / (4 * np.pi); H = -np.tensordot((2 * n + 1) / (n * (n + 1)) ** (m - 1), Ln, 1) / (4 * np.pi)
    Gi = np.linalg.inv(G + lam * np.eye(len(labs))); s = Gi.sum(0); C = Gi - np.outer(Gi.sum(1), s) / s.sum()
    return H @ C
SCALP = None; rows = []; Tc = {}; ST = {}
for stim in ("Real", "Sham"):
    for s in S17:
        for e in io.EYES:
            id_ = f"spTMS_{s:03d}_{stim}_{e}"
            if not (io.FINAL / f"{id_}_F_ERM.set").exists(): continue
            x, labs, tt, inten, _ = io.load_set(io.FINAL, id_ + "_F_ERM"); tt = np.round(tt).astype(int); w = (tt >= -1000) & (tt <= -10); assert w.sum() == 991
            if SCALP is None: SCALP = [l for l in labs if l not in ("M1", "M2")]
            key = tuple(labs)
            if key not in Tc: Tc[key] = csd_matrix(labs)
            seg = x[:, w, :].astype(float)
            for ver, dat in (("average reference", seg), ("Laplacian", np.einsum("ij,jtk->itk", Tc[key], seg))):
                f, Pw = signal.welch(dat, fs=1000, window="hann", nperseg=500, noverlap=250, axis=1); fa = (f >= 8) & (f <= 12)
                pa = Pw[:, fa, :].mean(1)                                   # channel x trial
                for c in SCALP: rows.append(dict(stim=stim, sub=s, eye=e, version=ver, channel=c, log10_power=np.log10(pa[labs.index(c)].mean())))
                if stim == "Real" and ver == "average reference": ST[id_] = (np.log10(pa[[labs.index(c) for c in SCALP]]).T.astype(np.float32), inten.astype(int))
            print(id_, flush=True)
A = pd.DataFrame(rows); A.to_csv(D / "EyeState_channel_alpha.csv", index=False)
T = []
for (stim, ver), g in A.groupby(["stim", "version"]):
    for c in SCALP:
        w_ = g[g.channel == c].pivot(index="sub", columns="eye", values="log10_power").dropna(); d = w_.EC - w_.EO; r = stats.ttest_rel(w_.EC, w_.EO)
        T.append(dict(stim=stim, version=ver, channel=c, n=len(d), t=r.statistic, p=r.pvalue, dz=d.mean() / d.std(ddof=1)))
T = pd.DataFrame(T)
from statsmodels.stats.multitest import multipletests
T["q_fdr61"] = np.nan
for k, g in T.groupby(["stim", "version"]): T.loc[g.index, "q_fdr61"] = multipletests(g.p, method="fdr_bh")[1]
T.to_csv(D / "EyeState_channel_alpha_t.csv", index=False)
ids = sorted(ST); np.savez_compressed(D / "Prestim_power_single_trial_allch.npz", ids=np.array(ids), channels=np.array(SCALP),
                                      **{f"P_{i}": ST[i][0] for i in ids}, **{f"I_{i}": ST[i][1] for i in ids})
# (3) grand-average conventional TEP
GA = {}
for e in io.EYES:
    subs = []
    for s in S17:
        x, labs, tt, inten, _ = io.load_set(io.FINAL, f"spTMS_{s:03d}_Real_{e}_F_TEP"); tt = np.round(tt).astype(int); w = (tt >= -100) & (tt <= 400)
        subs.append(np.stack([x[[labs.index(c) for c in SCALP]][:, w][:, :, inten == lv].mean(2) for lv in io.LEV]))
    GA[e] = np.stack(subs).mean(0)                                          # intensity x channel x time
np.savez_compressed(D / "TEP_GA_allch.npz", EO=GA["EO"], EC=GA["EC"], time=tt[w], channels=np.array(SCALP), levels=np.array(io.LEV))
print("done", T.groupby(["stim", "version"]).apply(lambda g: (g.q_fdr61 < .05).sum()).to_dict())
