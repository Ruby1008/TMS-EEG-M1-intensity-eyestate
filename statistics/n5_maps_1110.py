"""N5 per-channel maps (user 2026-10-28), as Figure 7C (code/round9_mu_stats.py, part D): for every scalp channel, single-trial log10 8-12 Hz power
(data/Prestim_power_single_trial_allch.npz) centred within participant x intensity; ROI P30 variant ~ power_within + C(intensity) + (1 | participant),
ML (statsmodels MixedLM); t of power_within; BH-FDR over 61 channels. Variants: main (must reproduce data/SingleTrial_mu_channel_P30.csv),
IMF4-5, IMF6, conventional TEP (from _check/N5_single_trial.csv). Figures: FigureS_muP30_IMF45_1110, FigureS_muP30_convTEP_1110 (layout of the
Figure 7C panel: coolwarm t map, FDR circles, sensorimotor ROI outline, colour bar 't'). Also writes N5/N5_mu_P30_variants.csv (merging _check/N5_lmm.csv)."""
data_dir = r"path/to/data"   # release: set this path
import warnings
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
from scipy.interpolate import griddata
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import r1110_lib as L
warnings.filterwarnings("ignore")
OUT = L.RV / "N5"
Z = np.load(L.DATA / "Prestim_power_single_trial_allch.npz"); CH = list(Z["channels"]); ids = list(Z["ids"])
ST = pd.read_csv(L.CK / "N5_single_trial.csv")
rows = []
for i in ids:
    P = Z[f"P_{i}"]; I = Z[f"I_{i}"]; st = ST[ST.id == i].sort_values("epoch"); assert len(st) == len(I) and np.array_equal(st.intensity.values, I)
    d = pd.DataFrame(P, columns=CH)
    for v in ("P30_main", "P30_IMF45", "P30_IMF6", "P30_convTEP"): d[v] = st[v].values
    d["intensity"] = I; d["sub"] = int(i.split("_")[1]); rows.append(d)
A = pd.concat(rows, ignore_index=True)
res = []
for v in ("P30_main", "P30_IMF45", "P30_IMF6", "P30_convTEP"):
    for c in CH:
        x = A[["sub", "intensity", v, c]].rename(columns={c: "pw", v: "y"}).copy(); x["pw_w"] = x.pw - x.groupby(["sub", "intensity"]).pw.transform("mean")
        m = smf.mixedlm("y ~ pw_w + C(intensity)", x, groups=x["sub"]).fit(reml=False)
        res.append(dict(variant=v, channel=c, b=m.params["pw_w"], se=m.bse["pw_w"], t=m.tvalues["pw_w"], p=m.pvalues["pw_w"]))
    print(v, flush=True)
M = pd.DataFrame(res)
for v, g in M.groupby("variant"): M.loc[g.index, "q_fdr61"] = L.bh(g.p)
M.to_csv(L.CK / "N5_channel_maps.csv", index=False)
REF = pd.read_csv(L.DATA / "SingleTrial_mu_channel_P30.csv").set_index("channel"); mm = M[M.variant == "P30_main"].set_index("channel").loc[REF.index]
assert np.allclose(mm.t, REF.t, atol=1e-4) and ((mm.q_fdr61 < .05) == (REF.q_fdr61 < .05)).all(), "main map does not reproduce Figure 7C"
NS = M.groupby("variant").apply(lambda g: int((g.q_fdr61 < .05).sum())).to_dict(); print("FDR channels:", NS)
# figures
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 9, "axes.linewidth": .8})
CED = pd.read_csv(data_dir + r"\waveguard64_woRef.ced", sep="\t", skipinitialspace=True); CED.columns = [c.strip() for c in CED.columns]; CED["labels"] = CED["labels"].str.strip()
POS = {r.labels: (float(r.radius) * np.sin(np.deg2rad(float(r.theta))), float(r.radius) * np.cos(np.deg2rad(float(r.theta)))) for r in CED.itertuples()}
def roi_outline(ax, grp, m=0.05):
    g = np.array([POS[l] for l in grp]); c = (g.max(0) + g.min(0)) / 2; w, h = g.max(0) - g.min(0) + 2 * m
    ax.add_patch(Ellipse(c, w * 1.25, h * 1.25 + 0.02, fill=False, ls="--", lw=.8, color="k"))
def topo(ax, labs, vals, vmax, sig=None):
    xy = np.array([POS[l] for l in labs]); R = np.hypot(*xy.T).max() * 1.03; ang = np.linspace(0, 2 * np.pi, 48, endpoint=False); ring = np.c_[R * np.cos(ang), R * np.sin(ang)]
    rv = np.array([vals[np.argsort(np.hypot(*(xy - r_).T))[:3]].mean() for r_ in ring]); gx, gy = np.meshgrid(np.linspace(-R, R, 300), np.linspace(-R, R, 300))
    z = griddata(np.r_[xy, ring], np.r_[vals, rv], (gx, gy), method="linear"); z[np.hypot(gx, gy) > R] = np.nan
    im = ax.imshow(z, extent=(-R, R, -R, R), origin="lower", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.plot(*POS["FC3"], marker=(6, 2, 0), color="k", ms=9, mew=1.6)
    if sig is not None:
        for l, s in zip(labs, sig):
            if s: ax.plot(*POS[l], "o", mfc="none", mec="k", ms=6.5, mew=1.1)
    ax.set_xlim(-R * 1.02, R * 1.02); ax.set_ylim(-R * 1.02, R * 1.02); ax.set_aspect("equal"); ax.axis("off"); return im
vm = np.ceil(M[M.variant.isin(["P30_main", "P30_IMF45", "P30_convTEP"])].t.abs().max())
for v, stem, ttl in (("P30_IMF45", "FigureS_muP30_IMF45_1110", "Channel mu power → ROI P30 (ERM, IMF4–5)"), ("P30_convTEP", "FigureS_muP30_convTEP_1110", "Channel mu power → ROI P30 (conventional TEP)")):
    g = M[M.variant == v]; fig = plt.figure(figsize=(3.0, 2.6)); ax = fig.add_axes([0.02, 0.04, 0.72, 0.82])
    im = topo(ax, list(g.channel), g.t.values, vm, sig=(g.q_fdr61 < .05).values); roi_outline(ax, list(L.io.ROI)); ax.set_title(ttl, fontsize=9, pad=8)
    cb = fig.add_axes([0.80, 0.18, 0.03, 0.55]); fig.colorbar(im, cax=cb); cb.set_title("t", fontsize=9); cb.tick_params(labelsize=7)
    for ext in ("png", "pdf"): fig.savefig(OUT / f"{stem}.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
# CSV
LM = pd.read_csv(L.CK / "N5_lmm.csv")
mapv = {"P30_main": "P30_main", "P30_IMF45": "P30_IMF45", "P30_IMF6": "P30_IMF6", "P30_convTEP": "P30_convTEP"}
LM["n_channels_sig"] = [NS.get(o) if p == "log_mu" else np.nan for o, p in zip(LM.outcome, LM.predictor)]
LM.rename(columns={"change_per_SD_uV": "change_per_SD_uV"})[["variant", "outcome", "predictor", "intensity", "n", "n_trials", "b", "SE", "df", "t", "p", "change_per_SD_uV", "n_channels_sig"]].to_csv(OUT / "N5_mu_P30_variants.csv", index=False)
print(M.groupby("variant").t.agg(["min", "max"]).round(2).to_string())
print(M[M.channel.isin(L.io.ROI)].pivot(index="channel", columns="variant", values="t").round(2).to_string())
