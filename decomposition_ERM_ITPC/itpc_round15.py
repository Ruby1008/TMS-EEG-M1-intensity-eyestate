"""Round 15 (user prompt 2026-10-27, section 1): Figure 5 (ITPC) with a new baseline and post-stimulus testing only.
New (main): baseline -500..-100 ms (subtract the mean ITPC over -500..-100 ms per channel-ROI x IMF x intensity; avoids the
  pre-pulse phase-estimation smearing zone); 2 (eye) x 3 (intensity) rmANOVA (GG) and post hoc paired t at IMF 3-8 x 0..500 ms
  only; BH-FDR q = .05 per effect (and per post hoc contrast) over IMF x 0..500 ms; -100..0 ms is plotted but not tested.
  No "value at 0 ms as zero point" correction.
Old (pre-specified, now Supplementary sensitivity): baseline -500..-10 ms, testing over -100..500 ms (= itpc_fig78_G.py; reproduced
  here and asserted identical to data/Fig7_ITPC_significant_segments.csv and Fig7_ITPC_posthoc_FDR_segments.csv).
Input (read only): data/ITPC_roi_real.mat (X = ROI ITPC baseline-corrected with -500..-10 ms; raw = X + constant, so
  re-baselining X with -500..-100 ms equals re-baselining the raw ITPC).
Statistics functions copied unchanged from itpc_fig78_G.py (rm2, bh).
Outputs: data/ITPC_r15_{new,old}_anova_maps.csv, ITPC_r15_{new,old}_significant_segments.csv, ITPC_r15_{new,old}_posthoc_segments.csv,
  ITPC_r15_baseline_comparison.csv (segment differences), ITPC_r15_longest_segments.csv;
  Fig5_r15.png/.tiff (main, new), FigS_ITPC_prespecified_r15.png/.tiff (old version, same layout).
Panel A reproduces the NeuroHolo2 perform_RmAnova2.m plotting calls (imagesc of F, white contour of the FDR mask, LineWidth 2,
  line at 0 ms, colorbar) with the colormap and 0-12 range of the original Figure 7 (perform_RmAnova2 needs an SPM MEEG object and
  computes its own uncorrected ANOVA, so it is not called directly)."""
from pathlib import Path
import numpy as np, pandas as pd, scipy.io as sio
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(__file__).resolve().parent.parent; D = ROOT / "data"
M = sio.loadmat(D / "ITPC_roi_real.mat", squeeze_me=True)
X = M["X"]; ids = [str(i) for i in M["ids"]]; tf = M["t"].astype(float); t = np.round(tf)   # tf: float axis as used by itpc_fig78_G.py
Q = pd.read_csv(ROOT.parent / "Preproc" / "QC" / "QC_summary_wide.csv", encoding="utf-8-sig")
main = Q[(Q.status == "processed") & (Q.stim == "Real") & (Q.stageD != "exclude")]
subs = sorted(set(main[main.eyes == "EO"]["sub"]) & set(main[main.eyes == "EC"]["sub"])); idx = {i: k for k, i in enumerate(ids)}
IMFS = list(range(3, 9)); tw = (t >= -100) & (t <= 500); tt = t[tw]                     # new: rounded axis, -100..500 (601 points)
twf = (tf >= -100) & (tf <= 500)   # old: float axis of itpc_fig78_G.py; -100 ms is -100.0000000001 there, so the old analysis used -99..500 (600 points)
XN = X - X[:, :, (t >= -500) & (t <= -100), :].mean(2, keepdims=True)          # new baseline -500..-100 ms

def stackY(XX, tw):
    Y = np.stack([np.stack([XX[idx[f"spTMS_{s:03d}_Real_{e}"]][np.array(IMFS) - 1][:, tw, :] for e in ("EO", "EC")], 0) for s in subs], 0)
    return np.moveaxis(Y, 4, 2)                                                    # n x eye x int x imf x time

def rm2(Y):
    n = Y.shape[0]
    GM = Y.mean((0, 1, 2)); mS = Y.mean((1, 2)); mA = Y.mean((0, 2)); mB = Y.mean((0, 1))
    mAS = Y.mean(2); mBS = Y.mean(1); mAB = Y.mean(0)
    a, b = 2, 3
    SSA = b * n * ((mA - GM) ** 2).sum(0); SSAS = b * ((mAS - mA[None] - mS[:, None] + GM) ** 2).sum((0, 1))
    SSB = a * n * ((mB - GM) ** 2).sum(0); SSBS = a * ((mBS - mB[None] - mS[:, None] + GM) ** 2).sum((0, 1))
    SSAB = n * ((mAB - mA[:, None] - mB[None] + GM) ** 2).sum((0, 1))
    SSABS = ((Y - mAB[None] - mAS[:, :, None] - mBS[:, None] + mA[None, :, None] + mB[None, None] + mS[:, None, None] - GM) ** 2).sum((0, 1, 2))
    out = {}
    for name, ss, sse, df in (("eye", SSA, SSAS, 1), ("intensity", SSB, SSBS, 2), ("interaction", SSAB, SSABS, 2)):
        dfe = df * (n - 1); F = (ss / df) / (sse / dfe)
        if df == 1: eps = np.ones_like(F)
        else:
            Z = mBS if name == "intensity" else (Y[:, 0] - Y[:, 1])
            Zc = Z - Z.mean(0)
            C = np.einsum("ni...,nj...->ij...", Zc, Zc) / (n - 1)
            Cd = C - C.mean(0, keepdims=True) - C.mean(1, keepdims=True) + C.mean((0, 1), keepdims=True)
            eps = np.trace(Cd) ** 2 / ((3 - 1) * (Cd ** 2).sum((0, 1)))
            eps = np.clip(eps, 1 / (3 - 1), 1)
        p = stats.f.sf(F, df * eps, dfe * eps)
        out[name] = dict(F=F, p=p, eps=eps, df=df, dfe=dfe)
    return out

def bh(p):
    pf = p.ravel(); o = np.argsort(pf); m = len(pf)
    q = np.empty(m); q[o] = np.minimum.accumulate((pf[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.minimum(q, 1).reshape(p.shape)

PAIRS = [(2, 0, "120%RMT vs. 80%RMT"), (2, 1, "120%RMT vs. 100%RMT"), (1, 0, "100%RMT vs. 80%RMT")]
def segs_of(mask):
    d = np.diff(np.r_[0, mask.astype(int), 0]); return list(zip(np.where(d == 1)[0], np.where(d == -1)[0] - 1))

def analyse(XX, twv, test):
    """twv: window over the 1501-point axis; test: boolean over the window (time points entering the FDR family); untested points get q = NaN, sig = False.
    Results are returned on the rounded -100..500 axis (tt); a point outside twv (old: -100 ms) is NaN / not significant."""
    Y = stackY(XX, twv); R = rm2(Y); T = test; ttv = np.round(tf[twv])
    for k in R:
        q = np.full_like(R[k]["p"], np.nan); q[:, T] = bh(R[k]["p"][:, T]); R[k]["q"] = q; R[k]["sig"] = np.nan_to_num(q, nan=1) < .05
    P = Y.mean(1); Tv, Pv = [], []
    for i, j, _ in PAIRS:
        r = stats.ttest_rel(P[:, i], P[:, j], axis=0); Tv.append(r.statistic); Pv.append(r.pvalue)
    Tv = np.array(Tv); Pv = np.array(Pv); phq = np.full_like(Pv, np.nan)
    for c in range(3): phq[c][:, T] = bh(Pv[c][:, T])
    ph = np.nan_to_num(phq, nan=1) < .05
    k0 = int(np.where(tt == ttv[0])[0][0]); assert np.array_equal(tt[k0:k0 + len(ttv)], ttv) and k0 + len(ttv) == len(tt)
    def pad(a, fill):                                                              # put the window onto the -100..500 axis
        if k0 == 0: return a
        return np.concatenate([np.full(a.shape[:-1] + (k0,), fill, dtype=a.dtype), a], -1)
    for k in R:
        for f_ in ("F", "p", "eps", "q"): R[k][f_] = pad(R[k][f_], np.nan)
        R[k]["sig"] = pad(R[k]["sig"], False)
    Tv = pad(Tv, np.nan); ph = pad(ph, False); Y = pad(Y, np.nan)
    seg = [(k, imf, tt[s0], tt[s1], s1 - s0 + 1, v["F"][a, s0:s1 + 1].max()) for k, v in R.items() for a, imf in enumerate(IMFS) for s0, s1 in segs_of(v["sig"][a])]
    phs = [(lab, imf, tt[s0], tt[s1], s1 - s0 + 1, Tv[c, a, s0:s1 + 1].mean()) for c, (_, _, lab) in enumerate(PAIRS) for a, imf in enumerate(IMFS) for s0, s1 in segs_of(ph[c, a])]
    maps = pd.DataFrame([(k, imf, tm, v["F"][a, b], v["eps"][a, b], v["p"][a, b], v["q"][a, b]) for k, v in R.items() for a, imf in enumerate(IMFS) for b, tm in enumerate(tt)],
                        columns=["effect", "IMF", "time_ms", "F", "epsGG", "pGG", "qBH"])
    return dict(Y=Y, R=R, T=Tv, ph=ph, seg=pd.DataFrame(seg, columns=["effect", "IMF", "start_ms", "end_ms", "n_points", "F_max"]),
                phs=pd.DataFrame(phs, columns=["contrast", "IMF", "start_ms", "end_ms", "n_points", "mean_t"]), maps=maps)

OLD = analyse(X, twf, np.ones(twf.sum(), bool)); NEW = analyse(XN, tw, tt >= 0)
# reproduction check of the pre-specified analysis against the files of itpc_fig78_G.py
for mine, ref in ((OLD["seg"], "Fig7_ITPC_significant_segments.csv"), (OLD["phs"], "Fig7_ITPC_posthoc_FDR_segments.csv")):
    r = pd.read_csv(D / ref); assert len(r) == len(mine) and np.allclose(r.iloc[:, 1:].values.astype(float), mine.iloc[:, 1:].values.astype(float)), ref
print("old version reproduces itpc_fig78_G.py outputs: OK")
for tag, V in (("new", NEW), ("old", OLD)):
    V["maps"].to_csv(D / f"ITPC_r15_{tag}_anova_maps.csv", index=False); V["seg"].to_csv(D / f"ITPC_r15_{tag}_significant_segments.csv", index=False)
    V["phs"].to_csv(D / f"ITPC_r15_{tag}_posthoc_segments.csv", index=False)

# ---- comparison: per effect / contrast x IMF, points in 0..500 and -100..-1, segments, points only in one version
post = tt >= 0; rows = []
def segtxt(mask): return "; ".join(f"{tt[a]:.0f}–{tt[b]:.0f}" for a, b in segs_of(mask)) or "none"
for kind, keys in (("ANOVA", ("eye", "intensity", "interaction")), ("post hoc", tuple(p[2] for p in PAIRS))):
    for key in keys:
        for a, imf in enumerate(IMFS):
            if kind == "ANOVA": mo, mn = OLD["R"][key]["sig"][a], NEW["R"][key]["sig"][a]
            else: c = [p[2] for p in PAIRS].index(key); mo, mn = OLD["ph"][c, a], NEW["ph"][c, a]
            rows.append(dict(analysis=kind, effect_or_contrast=key, IMF=imf,
                             old_segments_post_0_500=segtxt(mo & post), new_segments_post_0_500=segtxt(mn & post),
                             old_n_post=int((mo & post).sum()), new_n_post=int((mn & post).sum()),
                             old_n_pre_m100_m1=int((mo & ~post).sum()), new_n_pre_m100_m1="not tested",
                             only_old_post=int((mo & ~mn & post).sum()), only_new_post=int((mn & ~mo & post).sum()),
                             old_segments_pre=segtxt(mo & ~post)))
CMP = pd.DataFrame(rows); CMP.to_csv(D / "ITPC_r15_baseline_comparison.csv", index=False)

# longest post-stimulus segment (0..500 ms) per effect / contrast and IMF (the numbers quoted in the Results)
def longest(mask):
    s = segs_of(mask)
    if not s: return None
    a, b = max(s, key=lambda x: x[1] - x[0]); return tt[a], tt[b], b - a + 1
lr = []
for tag, V in (("new", NEW), ("old", OLD)):
    for key in ("intensity",):
        for a, imf in enumerate(IMFS):
            L = longest(V["R"][key]["sig"][a] & post); lr.append((tag, "ANOVA " + key, imf) + (L if L else (np.nan, np.nan, 0)))
    for c, (_, _, lab) in enumerate(PAIRS):
        for a, imf in enumerate(IMFS):
            L = longest(V["ph"][c, a] & post); lr.append((tag, lab, imf) + (L if L else (np.nan, np.nan, 0)))
LR = pd.DataFrame(lr, columns=["version", "test", "IMF", "start_ms", "end_ms", "n_points"]); LR.to_csv(D / "ITPC_r15_longest_segments.csv", index=False)

# ---- figures
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 9, "axes.linewidth": .8})
cF = LinearSegmentedColormap.from_list("wr", ["#ffffff", "#f6c99a", "#e8835a", "#c0392b", "#7a0000"])   # original Figure 7 (itpc_fig78_G.py)
def plabel(ax, s): ax.text(-0.16, 1.08, s, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
def outline(ax, mask):
    """white contour of the binary IMF x time mask along cell edges (perform_RmAnova2.m: contour(..., 'w', 'LineWidth', 2))"""
    if not mask.any(): return
    up = 20; m = np.repeat(mask.astype(float), up, 0)
    yy = np.linspace(IMFS[0] - .5, IMFS[-1] + .5, m.shape[0] + 1)[:-1] + .5 / up
    m = np.pad(m, 1); yy = np.r_[yy[0] - 1 / up, yy, yy[-1] + 1 / up]; xx = np.r_[tt[0] - 1, tt, tt[-1] + 1]
    ax.contour(xx, yy, m, levels=[.5], colors="white", linewidths=1.2)

def figure(V, stem, title, smear):
    fig = plt.figure(figsize=(7.2, 6.4)); fig.suptitle(title, fontsize=11, fontweight="bold", y=.995)
    ext = [tt[0] - .5, tt[-1] + .5, IMFS[-1] + .5, IMFS[0] - .5]
    for j, (k, ttl) in enumerate((("eye", "Eye state"), ("intensity", "TMS intensity"), ("interaction", "Interaction"))):
        ax = fig.add_axes([0.07 + j * 0.325, 0.60, 0.23, 0.30])
        im = ax.imshow(V["R"][k]["F"], aspect="auto", cmap=cF, vmin=0, vmax=12, extent=ext, interpolation="nearest")
        outline(ax, V["R"][k]["sig"])
        if smear:
            ax.axvspan(-100, 0, color="0.80", alpha=.75, lw=0, zorder=3)
            ax.text(-50, 5.5, "phase-estimation smearing", rotation=90, ha="center", va="center", fontsize=6.5, color="0.15", zorder=4)
        ax.axvline(0, color="k", lw=1.2, zorder=5); ax.set_title(ttl, fontsize=9); ax.set_xlabel("Time (ms)"); ax.set_ylabel("IMF" if j == 0 else "")
        ax.set_ylim(IMFS[-1] + .5, IMFS[0] - .5); ax.set_xlim(tt[0], tt[-1]); ax.set_yticks(IMFS); ax.set_xticks(range(-100, 501, 100)); ax.tick_params(axis="x", labelsize=7.5)
        cb = fig.colorbar(im, ax=ax, fraction=.08, pad=.03); cb.set_ticks(range(0, 13, 2)); cb.set_label("F", labelpad=1)
        if j == 0: plabel(ax, "A")
    hl = None
    for j, imf in enumerate((6, 7)):
        ax = fig.add_axes([0.07 + j * 0.49, 0.08, 0.41, 0.38]); top = 0; Y = V["Y"].mean(1)       # eye states pooled: n x int x imf x time
        for li, (lv, col) in enumerate(((80, "#9BB7B4"), (100, "#2E6E6A"), (120, "#C5550A"))):
            x = Y[:, li, IMFS.index(imf)]; m = x.mean(0); se = x.std(0, ddof=1) / np.sqrt(len(x))
            ax.fill_between(tt, m - se, m + se, color=col, alpha=.25, lw=0); ax.plot(tt, m, color=col, lw=1.4, label=f"{lv}% RMT"); top = max(top, np.nanmax(m + se))
        for k, (con, col) in enumerate((("120%RMT vs. 80%RMT", "#C5550A"), ("100%RMT vs. 80%RMT", "#2E6E6A"), ("120%RMT vs. 100%RMT", "0.3"))):
            yb = top * (1.06 + .06 * k); s_ = V["phs"][(V["phs"].contrast == con) & (V["phs"].IMF == imf)]
            if smear: s_ = s_[s_.start_ms >= 0]
            for r in s_.itertuples(): ax.plot([r.start_ms, r.end_ms], [yb, yb], color=col, lw=3, solid_capstyle="butt")
            ax.plot([], [], color=col, lw=3, label=con.replace("%RMT", "%").replace(" vs. ", " vs ") + ("" if len(s_) else " (none)"))
        ax.axvline(0, color="k", lw=.6, ls=":"); ax.axhline(0, color="0.6", lw=.5); ax.set_xlim(-100, 500); ax.set_xlabel("Time (ms)")
        ax.set_ylabel("ITPC (baseline-corrected)" if j == 0 else ""); ax.set_title(f"IMF{imf} ({'alpha' if imf == 6 else 'theta'})", fontsize=9)
        ax.text(0.03, 0.04, "80% < 100% ≈ 120%", transform=ax.transAxes, fontsize=9, fontweight="bold"); ax.set_ylim(top=top * 1.25)
        if j == 0: plabel(ax, "B"); hl = ax.get_legend_handles_labels()
    fig.legend(*hl, loc="lower center", bbox_to_anchor=(0.5, -0.06), ncol=6, frameon=False, fontsize=7.5)
    for ext_ in ("png", "tiff"):
        fig.savefig(ROOT / f"{stem}.{ext_}", dpi=600, bbox_inches="tight", **({"pil_kwargs": {"compression": "tiff_lzw"}} if ext_ == "tiff" else {}))
    plt.close(fig)

TITLE = "Phase locking rises from 80% to 100% RMT and does not increase further at 120% RMT"
figure(NEW, "Fig5_r15", TITLE, smear=True)
figure(OLD, "FigS_ITPC_prespecified_r15", "Pre-specified ITPC analysis (baseline −500 to −10 ms; tested from −100 to 500 ms)", smear=False)
pd.set_option("display.width", 250)
print(CMP.to_string()); print(LR.to_string())
