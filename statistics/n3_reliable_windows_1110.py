"""N3 (user 2026-10-28): is the 80% RMT reliable response a residual artefact?
_1110 copy of code/r21_round15.py (main analysis = 'local' null: each null segment minus the mean of the 94 ms before it; 2000 draws; seed 20261017;
single-trial RMS window vs -100..-55 ms, one-sided Wilcoxon; BH-FDR over all real and sham cells per test; responder = both q < .05).
Changes: the response window is a parameter. Windows run: 15-60 ms (main; reproduction check against data/R21_r15_cells.csv), 25-60, 30-60
and 15-25 ms (segment split for d). One RNG per window, all seeded 20261017 and consumed in the original cell order, so 15-60 reproduces exactly.
The null segment length equals the window length; the pre-stimulus RMS window stays -100..-55 ms.
(a) 80% RMT topographies (EO, EC): mean absolute amplitude 15-60 ms of the trial average per channel (real, sham, real - sham) and the paired
    t map real - sham (participants with both, BH-FDR over 61 channels). Head drawing: topo() of Figures_all_2026-11-08/code/fig6_a3_rev1108.py.
(c) residual (Preproc/QC/QC_summary_wide.csv residual11to25_<lv>_uV = max |ROI average| 11-25 ms, the stage-D metric) vs 15-60 ms log SNR (real):
    Spearman pooled (6 conditions x 17) and per intensity, repeated-measures correlation (pingouin.rm_corr, participants).
(d) 15-25 vs 25-60 ms: GMFP and log SNR, real vs sham.
Outputs (N3/): N3_reliable_windows.csv, N3_residual_vs_SNR.csv, N3_GMFP_segments.csv, FigureS_topo80_1110.png/.pdf; _check/N3_cells_all_windows.csv"""
data_dir = r"path/to/data"   # release: set this path
import numpy as np, pandas as pd, pingouin as pg
from scipy import stats
from scipy.interpolate import griddata
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import r1110_lib as L
io = L.io; OUT = L.RV / "N3"
B = 2000; LB = 94
WINS = {"15-60": (15, 60), "25-60": (25, 60), "30-60": (30, 60), "15-25": (15, 25)}
RNG = {w: np.random.default_rng(20261017) for w in WINS}
rows, topo80 = [], {}
for stim in ("Real", "Sham"):
    for s in io.MAIN17:
        for e in io.EYES:
            id_ = f"spTMS_{s:03d}_{stim}_{e}"
            if not (io.FINAL / f"{id_}_F_TEP.set").exists(): continue
            x, labs, t, inten, _ = io.load_set(io.FINAL, id_ + "_F_TEP"); t = np.round(t).astype(int)
            sc = [i for i, l in enumerate(labs) if l not in ("M1", "M2")]; x = x[sc].astype(float); slabs = [labs[i] for i in sc]
            starts = np.where((t >= -1000) & (t <= -56))[0]; pre_rms = np.where((t >= -100) & (t <= -55))[0]
            for lv in io.LEV:
                X = x[:, :, inten == lv]; n = X.shape[2]; Xt = np.transpose(X, (2, 1, 0))
                cs = np.concatenate([np.zeros((n, 1, Xt.shape[2])), np.cumsum(Xt, 1)], 1)
                rms_pre = np.sqrt((X[:, pre_rms, :] ** 2).mean((0, 1)))
                if lv == 80:
                    avg = X.mean(2); topo80[(stim, s, e)] = (slabs, np.abs(avg[:, (t >= 15) & (t <= 60)]).mean(1))
                for wn, (a0, a1) in WINS.items():
                    post = np.where((t >= a0) & (t <= a1))[0]; Lw = len(post); assert Lw == a1 - a0 + 1
                    obs = X[:, post, :].mean(2).std(0).mean(); null = np.empty(B); nloc = np.empty(B)
                    for b0 in range(0, B, 100):
                        st = RNG[wn].choice(starts, size=(100, n)); idx = st[:, :, None] + np.arange(Lw)
                        seg = Xt[np.arange(n)[None, :, None], idx]
                        null[b0:b0 + 100] = seg.mean(1).std(2).mean(1)
                        base = (cs[np.arange(n)[None, :], st] - cs[np.arange(n)[None, :], st - LB]) / LB
                        nloc[b0:b0 + 100] = (seg - base[:, :, None, :]).mean(1).std(2).mean(1)
                    rms_post = np.sqrt((X[:, post, :] ** 2).mean((0, 1)))
                    p_w = stats.wilcoxon(rms_post, rms_pre, alternative="greater").pvalue
                    rows.append(dict(window=wn, id=id_, sub=s, stim=stim, eye=e, intensity=lv, n_trials=n, gmfp_post=obs, p_rms=p_w,
                                     gmfp_null_mean=nloc.mean(), log_snr=np.log(obs / nloc.mean()), p_boot=(1 + (nloc >= obs).sum()) / (B + 1),
                                     log_snr_orignull=np.log(obs / null.mean())))
            print(id_, flush=True)
C = pd.DataFrame(rows)
for wn, g in C.groupby("window"):
    C.loc[g.index, "q_rms"] = L.bh(g.p_rms); C.loc[g.index, "q_boot"] = L.bh(g.p_boot)
C["responder"] = (C.q_boot < .05) & (C.q_rms < .05)
C.to_csv(L.CK / "N3_cells_all_windows.csv", index=False)
R0 = pd.read_csv(L.DATA / "R21_r15_cells.csv").set_index(["id", "intensity"]); c0 = C[C.window == "15-60"].set_index(["id", "intensity"]).loc[R0.index]
assert np.allclose(c0.log_snr, R0.log_snr_local) and np.allclose(c0.p_boot, R0.p_boot_local) and (c0.responder.values == R0.responder_local.values).all()
print("15-60 ms reproduces data/R21_r15_cells.csv (local null): OK")
# (b) table
out = []
for wn in WINS:
    g = C[C.window == wn]
    for e in io.EYES:
        for lv in io.LEV:
            a = g[(g.stim == "Real") & (g.eye == e) & (g.intensity == lv)].set_index("sub"); b = g[(g.stim == "Sham") & (g.eye == e) & (g.intensity == lv)].set_index("sub")
            sub = a.index.intersection(b.index); ra, rb = a.loc[sub, "responder"], b.loc[sub, "responder"]; n10, n01 = int((ra & ~rb).sum()), int((~ra & rb).sum())
            p_mc = stats.binomtest(n10, n10 + n01, .5).pvalue if n10 + n01 > 0 else 1.0
            w = stats.wilcoxon(a.loc[sub, "log_snr"], b.loc[sub, "log_snr"])
            for stim, d in (("real", a), ("sham", b)):
                out.append(dict(window=wn, stim=stim, eye=e, intensity=lv, n_reliable=int(d.responder.sum()), n_total=len(d), logSNR_mean=d.log_snr.mean(),
                                logSNR_SE=d.log_snr.std(ddof=1) / np.sqrt(len(d)), GMFP_mean_uV=d.gmfp_post.mean(), mcnemar_n_pairs=len(sub),
                                mcnemar_real_only=n10, mcnemar_sham_only=n01, mcnemar_p=p_mc, logSNR_wilcoxon_p=w.pvalue))
RW = pd.DataFrame(out)
for wn, g in RW.groupby("window"):
    idx = g[g.stim == "real"].index; RW.loc[idx, "logSNR_wilcoxon_p_holm"] = L.holm(g.loc[idx, "logSNR_wilcoxon_p"]); RW.loc[g[g.stim == "sham"].index, "logSNR_wilcoxon_p_holm"] = RW.loc[idx, "logSNR_wilcoxon_p_holm"].values
RW.to_csv(OUT / "N3_reliable_windows.csv", index=False)
# (d) segments
seg = RW[RW.window.isin(["15-25", "25-60", "15-60"])][["window", "stim", "eye", "intensity", "n_total", "GMFP_mean_uV", "logSNR_mean", "logSNR_SE", "n_reliable"]].copy()
dd = []
for wn in ("15-25", "25-60", "15-60"):
    g = C[C.window == wn]
    for e in io.EYES:
        for lv in io.LEV:
            a = g[(g.stim == "Real") & (g.eye == e) & (g.intensity == lv)].set_index("sub"); b = g[(g.stim == "Sham") & (g.eye == e) & (g.intensity == lv)].set_index("sub")
            sub = a.index.intersection(b.index)
            dd.append(dict(window=wn, eye=e, intensity=lv, n_pairs=len(sub), GMFP_real_minus_sham=(a.loc[sub, "gmfp_post"] - b.loc[sub, "gmfp_post"]).mean(),
                           GMFP_ratio_real_over_sham=(a.loc[sub, "gmfp_post"] / b.loc[sub, "gmfp_post"]).median(),
                           logSNR_real_minus_sham=(a.loc[sub, "log_snr"] - b.loc[sub, "log_snr"]).mean(), p_wilcoxon_logSNR=stats.wilcoxon(a.loc[sub, "log_snr"], b.loc[sub, "log_snr"]).pvalue,
                           p_wilcoxon_GMFP=stats.wilcoxon(a.loc[sub, "gmfp_post"], b.loc[sub, "gmfp_post"]).pvalue))
DD = pd.DataFrame(dd); seg.to_csv(OUT / "N3_GMFP_segments.csv", index=False); DD.to_csv(OUT / "N3_GMFP_segments_real_minus_sham.csv", index=False)
# (c) residual vs log SNR
Q = io.Q[(io.Q.stim == "Real") & (io.Q.status == "processed")].set_index("id")
g = C[(C.window == "15-60") & (C.stim == "Real")].copy(); g["residual_uV"] = [Q.loc[i, f"residual11to25_{lv}_uV"] for i, lv in zip(g.id, g.intensity)]
rr = []
r = stats.spearmanr(g.residual_uV, g.log_snr); rr.append(dict(analysis="Spearman pooled (6 conditions)", n=len(g), r=r.statistic, p=r.pvalue))
for lv in io.LEV:
    h = g[g.intensity == lv]; r = stats.spearmanr(h.residual_uV, h.log_snr); rr.append(dict(analysis=f"Spearman {lv}% RMT (EO + EC)", n=len(h), r=r.statistic, p=r.pvalue))
rm = pg.rm_corr(data=g, x="residual_uV", y="log_snr", subject="sub"); rr.append(dict(analysis="repeated-measures correlation (participants)", n=len(g), r=rm["r"].iloc[0], p=rm["pval"].iloc[0], df=rm["dof"].iloc[0]))
for lv in io.LEV:
    h = g[g.intensity == lv]; rm = pg.rm_corr(data=h, x="residual_uV", y="log_snr", subject="sub"); rr.append(dict(analysis=f"repeated-measures correlation {lv}% RMT (EO, EC)", n=len(h), r=rm["r"].iloc[0], p=rm["pval"].iloc[0], df=rm["dof"].iloc[0]))
RR = pd.DataFrame(rr); RR.to_csv(OUT / "N3_residual_vs_SNR.csv", index=False)
g[["sub", "eye", "intensity", "residual_uV", "log_snr", "gmfp_post"]].to_csv(L.CK / "N3_residual_vs_SNR_values.csv", index=False)
desc_res = g.groupby("intensity").residual_uV.describe()

# (a) topographies
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 9, "axes.linewidth": .8})
CED = pd.read_csv(data_dir + r"\waveguard64_woRef.ced", sep="\t", skipinitialspace=True); CED.columns = [c.strip() for c in CED.columns]; CED["labels"] = CED["labels"].str.strip()
POS = {r.labels: (float(r.radius) * np.sin(np.deg2rad(float(r.theta))), float(r.radius) * np.cos(np.deg2rad(float(r.theta)))) for r in CED.itertuples()}
def topo(ax, labs, vals, vmin, vmax, cmap, sig=None):
    xy = np.array([POS[l] for l in labs]); R = np.hypot(*xy.T).max() * 1.03; ang = np.linspace(0, 2 * np.pi, 48, endpoint=False); ring = np.c_[R * np.cos(ang), R * np.sin(ang)]
    rv = np.array([vals[np.argsort(np.hypot(*(xy - r_).T))[:3]].mean() for r_ in ring]); gx, gy = np.meshgrid(np.linspace(-R, R, 300), np.linspace(-R, R, 300))
    z = griddata(np.r_[xy, ring], np.r_[vals, rv], (gx, gy), method="linear"); z[np.hypot(gx, gy) > R] = np.nan
    im = ax.imshow(z, extent=(-R, R, -R, R), origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.plot(*np.c_[R * np.cos(np.linspace(0, 2 * np.pi, 200)), R * np.sin(np.linspace(0, 2 * np.pi, 200))].T, color="k", lw=.6)
    ax.plot(*POS["FC3"], marker=(6, 2, 0), color="k", ms=7, mew=1.3)
    if sig is not None:
        for l, s_ in zip(labs, sig):
            if s_: ax.plot(*POS[l], "o", mfc="none", mec="k", ms=4.5, mew=.9)
    ax.set_xlim(-R * 1.02, R * 1.02); ax.set_ylim(-R * 1.02, R * 1.02); ax.set_aspect("equal"); ax.axis("off"); return im
fig, axs = plt.subplots(2, 4, figsize=(7.4, 4.0)); TOPO = []
amax = max(np.max(v[1]) for v in topo80.values())
for i, e in enumerate(io.EYES):
    sub = [s for s in io.MAIN17 if ("Real", s, e) in topo80 and ("Sham", s, e) in topo80]
    labs = topo80[("Real", sub[0], e)][0]; assert all(topo80[(st, s, e)][0] == labs for st in ("Real", "Sham") for s in sub)
    R_ = np.array([topo80[("Real", s, e)][1] for s in sub]); S_ = np.array([topo80[("Sham", s, e)][1] for s in sub]); D_ = R_ - S_
    tt_ = stats.ttest_rel(R_, S_, axis=0); q = L.bh(tt_.pvalue)
    for k, l in enumerate(labs): TOPO.append(dict(eye=e, channel=l, n=len(sub), real_mean_abs_uV=R_[:, k].mean(), sham_mean_abs_uV=S_[:, k].mean(), diff_uV=D_[:, k].mean(), t=tt_.statistic[k], p=tt_.pvalue[k], q_fdr=q[k]))
    vm = np.ceil(max(R_.mean(0).max(), S_.mean(0).max()))
    im1 = topo(axs[i, 0], labs, R_.mean(0), 0, vm, "Reds"); topo(axs[i, 1], labs, S_.mean(0), 0, vm, "Reds")
    dm = np.ceil(np.abs(D_.mean(0)).max()); im2 = topo(axs[i, 2], labs, D_.mean(0), -dm, dm, "coolwarm")
    tm = np.ceil(np.abs(tt_.statistic).max()); im3 = topo(axs[i, 3], labs, tt_.statistic, -tm, tm, "coolwarm", sig=q < .05)
    for ax_, im_, lab in ((axs[i, 1], im1, "Mean |amplitude| (µV)"), (axs[i, 2], im2, "Difference (µV)"), (axs[i, 3], im3, "t (real − sham)")):
        cb = fig.colorbar(im_, ax=ax_, fraction=.046, pad=.02, shrink=.8); cb.set_label(lab, fontsize=7); cb.ax.tick_params(labelsize=7)
    axs[i, 0].text(-0.25, 0.5, {"EO": "Eyes open", "EC": "Eyes closed"}[e] + "\n80% RMT", transform=axs[i, 0].transAxes, rotation=90, ha="center", va="center", fontsize=9)
for j, ttl in enumerate(("Real", "Sham", "Real − sham", "Paired t")): axs[0, j].set_title(ttl, fontsize=9)
fig.suptitle("Mean absolute amplitude 15–60 ms at 80% RMT", fontsize=10, fontweight="bold", y=1.0)
fig.tight_layout()
for ext in ("png", "pdf"): fig.savefig(OUT / f"FigureS_topo80_1110.{ext}", dpi=600, bbox_inches="tight")
plt.close(fig)
TP = pd.DataFrame(TOPO); TP.to_csv(OUT / "N3_topo80_channels.csv", index=False)
pd.set_option("display.width", 250)
print(RW.round(4).to_string()); print(DD.round(4).to_string()); print(RR.round(4).to_string()); print(desc_res.round(2).to_string())
print(TP.groupby("eye").apply(lambda h: pd.Series(dict(n_sig=(h.q_fdr < .05).sum(), max_t=h.t.max(), min_t=h.t.min(), top=", ".join(h.sort_values("t", ascending=False).channel.head(6))))).to_string())
