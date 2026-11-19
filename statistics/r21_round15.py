"""Round 15 (user prompt 2026-10-27, section 2.2): diagnosis of the negative sham log SNR in Figure 4B and a local-baseline null.
Same data, recordings, seed and random draws as r21_sham.py (main analysis, unchanged); for every bootstrap draw two null values are computed:
  orig : as r21_sham.py (null segments taken from the -100..-6 ms baseline-corrected epochs; asserted identical to data/R21_cells.csv)
  local: each trial's null segment minus that trial's mean over the 94 ms immediately before the segment (per channel), i.e. the
         segments get a local baseline like the observed 15-60 ms window, which lies right after the -100..-6 ms baseline.
The observed 15-60 ms GMFP, the RMS test, BH-FDR over all cells and the responder rule (both q < .05) are unchanged.
Diagnosis (per cell): mean single-trial RMS over channels in 100-ms bins from -1000 to -100 ms, in -100..-55 and in 15..60 ms
(distance from the baseline window), and the null mean of both versions.
Outputs: data/R21_r15_cells.csv (both versions), data/R21_r15_group_tests.csv (local version, same tests as r21_sham.py),
         data/R21_r15_diagnosis.csv, data/R21_r15_logsnr_summary.csv"""
import numpy as np, pandas as pd, pingouin as pg
from scipy import stats
import tep_io as io

RNG = np.random.default_rng(20261017); B = 2000; L = 46; LB = 94
rows, diag = [], []
for stim in ("Real", "Sham"):
    for s in io.MAIN17:
        for e in io.EYES:
            id_ = f"spTMS_{s:03d}_{stim}_{e}"
            if not (io.FINAL / f"{id_}_F_TEP.set").exists(): continue
            x, labs, t, inten, _ = io.load_set(io.FINAL, id_ + "_F_TEP"); t = np.round(t).astype(int)
            sc = [i for i, l in enumerate(labs) if l not in ("M1", "M2")]; x = x[sc].astype(float)
            post = np.where((t >= 15) & (t <= 60))[0]; assert len(post) == L
            starts = np.where((t >= -1000) & (t <= -56))[0]; pre_rms = np.where((t >= -100) & (t <= -55))[0]
            assert starts.min() - LB >= 0
            for lv in io.LEV:
                X = x[:, :, inten == lv]; n = X.shape[2]
                obs = X[:, post, :].mean(2).std(0).mean()
                Xt = np.transpose(X, (2, 1, 0))                                 # trial x time x chan
                cs = np.concatenate([np.zeros((n, 1, Xt.shape[2])), np.cumsum(Xt, 1)], 1)   # cs[:, k] = sum of samples < k
                null, nloc = np.empty(B), np.empty(B)
                for b0 in range(0, B, 100):
                    st = RNG.choice(starts, size=(100, n))
                    idx = st[:, :, None] + np.arange(L)
                    seg = Xt[np.arange(n)[None, :, None], idx]                   # draws x trials x L x chan
                    null[b0:b0 + 100] = seg.mean(1).std(2).mean(1)
                    base = (cs[np.arange(n)[None, :], st] - cs[np.arange(n)[None, :], st - LB]) / LB   # draws x trials x chan: mean of samples st-94..st-1
                    nloc[b0:b0 + 100] = (seg - base[:, :, None, :]).mean(1).std(2).mean(1)
                rms_post = np.sqrt((X[:, post, :] ** 2).mean((0, 1))); rms_pre = np.sqrt((X[:, pre_rms, :] ** 2).mean((0, 1)))
                p_w = stats.wilcoxon(rms_post, rms_pre, alternative="greater").pvalue
                rows.append(dict(id=id_, sub=s, stim=stim, eye=e, intensity=lv, n_trials=n, gmfp_post=obs, p_rms=p_w,
                                 gmfp_null_mean_orig=null.mean(), log_snr_orig=np.log(obs / null.mean()), p_boot_orig=(1 + (null >= obs).sum()) / (B + 1),
                                 gmfp_null_mean_local=nloc.mean(), log_snr_local=np.log(obs / nloc.mean()), p_boot_local=(1 + (nloc >= obs).sum()) / (B + 1)))
                d = dict(id=id_, sub=s, stim=stim, eye=e, intensity=lv)
                st_rms = np.sqrt((X ** 2).mean(0))                               # time x trial: single-trial RMS over channels
                for a in range(-1000, -100, 100): d[f"rms_{a}_{a + 99}"] = st_rms[(t >= a) & (t <= a + 99)].mean()
                d["rms_-100_-55"] = st_rms[pre_rms].mean(); d["rms_15_60"] = st_rms[post].mean()
                d["gmfp_avg_-100_-55"] = X[:, pre_rms, :].mean(2).std(0).mean()  # GMFP of the trial average inside the baseline window
                diag.append(d)
            print(id_, flush=True)
C = pd.DataFrame(rows)
def bh(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); q = np.empty(m); q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]; return np.minimum(q, 1)
C["q_rms"] = bh(C.p_rms)
for v in ("orig", "local"):
    C[f"q_boot_{v}"] = bh(C[f"p_boot_{v}"]); C[f"responder_{v}"] = (C[f"q_boot_{v}"] < .05) & (C.q_rms < .05)
# reproduction check against the main analysis
R0 = pd.read_csv(io.DATA / "R21_cells.csv").set_index(["id", "intensity"]); c0 = C.set_index(["id", "intensity"]).loc[R0.index]
assert np.allclose(c0.log_snr_orig, R0.log_snr) and np.allclose(c0.p_boot_orig, R0.p_boot) and (c0.responder_orig.values == R0.responder.values).all()
print("orig version reproduces data/R21_cells.csv: OK")
C.to_csv(io.DATA / "R21_r15_cells.csv", index=False); pd.DataFrame(diag).to_csv(io.DATA / "R21_r15_diagnosis.csv", index=False)

def group_tests(C, v):
    G = []; rsp, lsn = f"responder_{v}", f"log_snr_{v}"
    for (stim, e, lv), g in C.groupby(["stim", "eye", "intensity"]): G.append(dict(version=v, test="responders", stim=stim, eye=e, intensity=lv, n=len(g), k=int(g[rsp].sum()), value=g[rsp].mean()))
    for (stim, e), g in C.groupby(["stim", "eye"]):
        W = g.pivot(index="sub", columns="intensity", values=rsp).dropna().astype(int)
        k = W.sum(1); cj = W.sum(0); N = W.values.sum(); kk = 3
        Q = (kk - 1) * (kk * (cj ** 2).sum() - N ** 2) / (kk * N - (k ** 2).sum()) if (kk * N - (k ** 2).sum()) > 0 else np.nan
        G.append(dict(version=v, test="Cochran Q (intensity)", stim=stim, eye=e, n=len(W), value=Q, p=stats.chi2.sf(Q, kk - 1) if np.isfinite(Q) else np.nan))
    mc = []
    for e in io.EYES:
        for lv in io.LEV:
            a = C[(C.stim == "Real") & (C.eye == e) & (C.intensity == lv)].set_index("sub"); b = C[(C.stim == "Sham") & (C.eye == e) & (C.intensity == lv)].set_index("sub")
            sub = a.index.intersection(b.index); ra, rb = a.loc[sub, rsp], b.loc[sub, rsp]
            n10, n01 = int((ra & ~rb).sum()), int((~ra & rb).sum())
            p_mc = stats.binomtest(n10, n10 + n01, .5).pvalue if n10 + n01 > 0 else 1.0
            w = stats.wilcoxon(a.loc[sub, lsn], b.loc[sub, lsn])
            mc.append(dict(eye=e, intensity=lv, n_pairs=len(sub), real_only=n10, sham_only=n01, p_mcnemar=p_mc, logsnr_real=a.loc[sub, lsn].mean(), logsnr_sham=b.loc[sub, lsn].mean(), W=w.statistic, p_logsnr=w.pvalue))
    MC = pd.DataFrame(mc); MC["p_logsnr_holm"] = pg.multicomp(MC.p_logsnr.values, method="holm")[1]
    for r in MC.itertuples():
        G.append(dict(version=v, test="McNemar real vs sham", eye=r.eye, intensity=r.intensity, n=r.n_pairs, value=f"{r.real_only}/{r.sham_only}", p=r.p_mcnemar))
        G.append(dict(version=v, test="log SNR real vs sham (Wilcoxon)", eye=r.eye, intensity=r.intensity, n=r.n_pairs, value=f"{r.logsnr_real:.3f} vs {r.logsnr_sham:.3f}", W=r.W, p=r.p_logsnr, p_holm=r.p_logsnr_holm))
    return pd.DataFrame(G)
GT = pd.concat([group_tests(C, "orig"), group_tests(C, "local")]); GT.to_csv(io.DATA / "R21_r15_group_tests.csv", index=False)
S = C.melt(id_vars=["stim", "eye", "intensity"], value_vars=["log_snr_orig", "log_snr_local"], var_name="version").groupby(["version", "stim", "eye", "intensity"]).value.agg(["mean", "sem", "count"]).reset_index()
S.to_csv(io.DATA / "R21_r15_logsnr_summary.csv", index=False)
pd.set_option("display.width", 220); print(S.to_string()); print(GT.to_string())
print(pd.DataFrame(diag).groupby(["stim", "intensity"]).mean(numeric_only=True).drop(columns="sub").round(3).T.to_string())
