"""Round 15 merge (user 2026-10-22), unified ITPC rule 'v2': baseline -500..-100 ms, IMF3-7 (IMF8 excluded: intensity differences inside the
baseline window, i.e. non-stationary baseline), tests at 10..500 ms (end of the interpolated window), BH-FDR q = .05 with one family per effect
(and per post hoc contrast) over IMF3-7 x 10..500 ms. Statistics functions (rm2, bh) from itpc_round15.py (= itpc_fig78_G.py).
Verification of the pre-stimulus interval: itpc_r15_prestim_check_v2.py (0 points for all three effects).
Outputs (data/): ITPC_v2_{real,aepr}_anova_maps.csv, _significant_segments.csv, _posthoc_segments.csv, ITPC_v2_longest_segments.csv,
  ITPC_v2_windows.csv (longest intensity segment >= 20 ms, IMF6/7), ITPC_v2_jaccard.csv (real vs AEPr intensity masks, IMF6/7, 10..500 ms),
  ITPC_v2_vs_prespecified.csv (longest segments, v2 vs the pre-specified analysis), Fig8_ITPC_N100_input_v2.csv (17, longest segment),
  Fig8_ITPC_N100_input_allsigpoints_v2.csv, Fig8_ITPC_N100_input_allavailable_v2.csv, Fig8_ITPC_N100_input_10300_v2.csv (10-300 ms);
  figures Fig5_v2.png/.tiff, FigS22_ITPC_AEPr_v2.png/.tiff (layout of Fig5_v2)."""
from pathlib import Path
import numpy as np, pandas as pd, scipy.io as sio
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

CODE = Path(__file__).resolve().parent; ROOT = CODE.parent; D = ROOT / "data"
_src = (CODE / "itpc_round15.py").read_text(encoding="utf-8"); ns = {"__file__": str(CODE / "itpc_round15.py")}; exec(_src.split("OLD = analyse(")[0], ns)
rm2, bh, segs_of, subs, main = ns["rm2"], ns["bh"], ns["segs_of"], ns["subs"], ns["main"]
IMFS = [3, 4, 5, 6, 7]; T0 = 10
PAIRS = [(2, 0, "120%RMT vs. 80%RMT"), (2, 1, "120%RMT vs. 100%RMT"), (1, 0, "100%RMT vs. 80%RMT")]

def load(f):
    M = sio.loadmat(D / f, squeeze_me=True); X = M["X"]; ids = [str(i) for i in M["ids"]]; t = np.round(M["t"].astype(float))
    return X - X[:, :, (t >= -500) & (t <= -100), :].mean(2, keepdims=True), {i: k for k, i in enumerate(ids)}, t
def analyse(XN, idx, t):
    tw = (t >= -100) & (t <= 500); tt = t[tw]; test = tt >= T0
    Y = np.stack([np.stack([XN[idx[f"spTMS_{s:03d}_Real_{e}"]][np.array(IMFS) - 1][:, tw, :] for e in ("EO", "EC")], 0) for s in subs], 0); Y = np.moveaxis(Y, 4, 2)
    R = rm2(Y)
    for k in R:
        q = np.full_like(R[k]["p"], np.nan); q[:, test] = bh(R[k]["p"][:, test]); R[k]["q"] = q; R[k]["sig"] = np.nan_to_num(q, nan=1) < .05
    P = Y.mean(1); Tv, Pv = [], []
    for i, j, _ in PAIRS:
        r = stats.ttest_rel(P[:, i], P[:, j], axis=0); Tv.append(r.statistic); Pv.append(r.pvalue)
    Tv, Pv = np.array(Tv), np.array(Pv); ph = np.zeros_like(Pv, bool)
    for c in range(3): q = bh(Pv[c][:, test]); ph[c][:, test] = q < .05
    seg = pd.DataFrame([(k, imf, tt[a], tt[b], b - a + 1, v["F"][i, a:b + 1].max()) for k, v in R.items() for i, imf in enumerate(IMFS) for a, b in segs_of(v["sig"][i])],
                       columns=["effect", "IMF", "start_ms", "end_ms", "n_points", "F_max"])
    phs = pd.DataFrame([(lab, imf, tt[a], tt[b], b - a + 1, Tv[c, i, a:b + 1].mean()) for c, (_, _, lab) in enumerate(PAIRS) for i, imf in enumerate(IMFS) for a, b in segs_of(ph[c, i])],
                       columns=["contrast", "IMF", "start_ms", "end_ms", "n_points", "mean_t"])
    maps = pd.DataFrame([(k, imf, tm, v["F"][i, b], v["eps"][i, b], v["p"][i, b], v["q"][i, b]) for k, v in R.items() for i, imf in enumerate(IMFS) for b, tm in enumerate(tt)],
                        columns=["effect", "IMF", "time_ms", "F", "epsGG", "pGG", "qBH"])
    return dict(Y=Y, R=R, ph=ph, tt=tt, seg=seg, phs=phs, maps=maps)
def longest(mask, tt):
    s = segs_of(mask)
    if not s: return None
    a, b = max(s, key=lambda x: x[1] - x[0]); return tt[a], tt[b], b - a + 1

XR, idxR, t = load("ITPC_roi_real.mat"); XA, idxA, tA = load("ITPC_roi_real_AEPr.mat"); assert np.array_equal(t, tA)
V = {"real": analyse(XR, idxR, t), "aepr": analyse(XA, idxA, t)}
for tag, v in V.items():
    for k in ("maps", "seg", "phs"): v[k].to_csv(D / f"ITPC_v2_{tag}_{ {'maps': 'anova_maps', 'seg': 'significant_segments', 'phs': 'posthoc_segments'}[k]}.csv".replace("_ ", "_").replace(" .", "."), index=False)
    for e in ("eye", "interaction"): assert not v["R"][e]["sig"].any() or tag == "aepr", (tag, e, int(v["R"][e]["sig"].sum()))
lr = []
for tag, v in V.items():
    tt = v["tt"]
    for i, imf in enumerate(IMFS):
        L = longest(v["R"]["intensity"]["sig"][i], tt); lr.append((tag, "ANOVA intensity", imf) + (L if L else (np.nan, np.nan, 0)))
        for e in ("eye", "interaction"):
            L = longest(v["R"][e]["sig"][i], tt); lr.append((tag, f"ANOVA {e}", imf) + (L if L else (np.nan, np.nan, 0)))
        for c, (_, _, lab) in enumerate(PAIRS):
            L = longest(v["ph"][c, i], tt); lr.append((tag, lab, imf) + (L if L else (np.nan, np.nan, 0)))
LR = pd.DataFrame(lr, columns=["data", "test", "IMF", "start_ms", "end_ms", "n_points"]); LR.to_csv(D / "ITPC_v2_longest_segments.csv", index=False)
# windows for the summary measures (rule of itpc_fig78_G.py: longest continuous segment >= 20 ms; all significant points for the S-H version)
WIN, wrow = {}, []
tt = V["real"]["tt"]
for imf in (6, 7):
    sig = V["real"]["R"]["intensity"]["sig"][IMFS.index(imf)]; L = longest(sig, tt); m_ = np.zeros_like(sig)
    assert L is not None and L[2] >= 20
    a = int(np.where(tt == L[0])[0][0]); m_[a:a + L[2]] = True; WIN[imf] = dict(main=m_, allp=sig.copy())
    wrow.append((imf, L[0], L[1], L[2], int(sig.sum()), f"{tt[sig].min():.0f}-{tt[sig].max():.0f}"))
pd.DataFrame(wrow, columns=["IMF", "longest_start_ms", "longest_end_ms", "longest_n_points", "all_sig_n_points", "all_sig_range_ms"]).to_csv(D / "ITPC_v2_windows.csv", index=False)
# Jaccard of the intensity masks, real vs AEPr (IMF6, IMF7; 10..500 ms)
J = []
for imf in (6, 7):
    a, b = V["real"]["R"]["intensity"]["sig"][IMFS.index(imf)], V["aepr"]["R"]["intensity"]["sig"][IMFS.index(imf)]
    J.append(dict(IMF=imf, n_real=int(a.sum()), n_aepr=int(b.sum()), n_both=int((a & b).sum()), jaccard=(a & b).sum() / (a | b).sum()))
pd.DataFrame(J).to_csv(D / "ITPC_v2_jaccard.csv", index=False)
# comparison with the pre-specified analysis (data/ITPC_r15_longest_segments.csv, version 'old')
OLD = pd.read_csv(D / "ITPC_r15_longest_segments.csv"); OLD = OLD[OLD.version == "old"]
cmp_ = []
for r in LR[(LR.data == "real") & ~LR.test.isin(["ANOVA eye", "ANOVA interaction"])].itertuples():
    o = OLD[(OLD.test == r.test) & (OLD.IMF == r.IMF)].iloc[0]
    f = lambda a, b, n: f"{a:.0f}–{b:.0f}" if n else "none"
    cmp_.append(dict(test=r.test, IMF=r.IMF, prespecified=f(o.start_ms, o.end_ms, o.n_points), v2=f(r.start_ms, r.end_ms, r.n_points)))
for imf in (8,):
    for test in ("ANOVA intensity", "120%RMT vs. 80%RMT", "120%RMT vs. 100%RMT", "100%RMT vs. 80%RMT"):
        o = OLD[(OLD.test == test) & (OLD.IMF == imf)].iloc[0]; cmp_.append(dict(test=test, IMF=8, prespecified=f"{o.start_ms:.0f}–{o.end_ms:.0f}" if o.n_points else "none", v2="not analyzed"))
pd.DataFrame(cmp_).to_csv(D / "ITPC_v2_vs_prespecified.csv", index=False)

# ---- per-participant ITPC window means (as itpc_fig78_G.py f8table / fig8_window0300.py), N100 = ERM ROI 120-140 ms
HH = ROOT / "HHSA"; tw = (t >= -100) & (t <= 500)
def f8table(subjects, kind=None, fixed=None):
    rows = []
    for s_ in subjects:
        r = {"sub": s_}
        for eye in ("EO", "EC"):
            id_ = f"spTMS_{s_:03d}_Real_{eye}"; ok = id_ in idxR and ((main["sub"] == s_) & (main.eyes == eye)).any()
            if ok: E = sio.loadmat(HH / f"erm_{id_}.mat", squeeze_me=True); te = np.round(E["time"])
            for k, lv in enumerate((80, 100, 120)):
                for imf in (6, 7):
                    w = WIN[imf][kind] if kind else (tt >= fixed[0]) & (tt <= fixed[1])
                    r[f"ITPC_IMF{imf}_{eye}_{lv}"] = XR[idxR[id_], imf - 1][tw][w, k].mean() if ok else np.nan
                r[f"N100_{eye}_{lv}"] = E["erm_roi"][k][(te >= 120) & (te <= 140)].mean() if ok else np.nan
        rows.append(r)
    return pd.DataFrame(rows)
f8table(subs, "main").to_csv(D / "Fig8_ITPC_N100_input_v2.csv", index=False)
f8table(subs, "allp").to_csv(D / "Fig8_ITPC_N100_input_allsigpoints_v2.csv", index=False)
f8table(sorted(set(main["sub"])), "main").to_csv(D / "Fig8_ITPC_N100_input_allavailable_v2.csv", index=False)
f8table(subs, fixed=(10, 300)).to_csv(D / "Fig8_ITPC_N100_input_10300_v2.csv", index=False)
# check: N100 values identical to the existing input
o8 = pd.read_csv(D / "Fig8_ITPC_N100_input.csv"); n8 = pd.read_csv(D / "Fig8_ITPC_N100_input_v2.csv")
assert np.allclose(o8.filter(like="N100").values, n8.filter(like="N100").values)

# ---- figures (layout of Fig5_r15; y axis IMF3-7)
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 9, "axes.linewidth": .8})
cF = LinearSegmentedColormap.from_list("wr", ["#ffffff", "#f6c99a", "#e8835a", "#c0392b", "#7a0000"])
def plabel(ax, s): ax.text(-0.16, 1.08, s, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
def outline(ax, mask, tt):
    if not mask.any(): return
    up = 20; m = np.repeat(mask.astype(float), up, 0)
    yy = np.linspace(IMFS[0] - .5, IMFS[-1] + .5, m.shape[0] + 1)[:-1] + .5 / up
    m = np.pad(m, 1); yy = np.r_[yy[0] - 1 / up, yy, yy[-1] + 1 / up]; xx = np.r_[tt[0] - 1, tt, tt[-1] + 1]
    ax.contour(xx, yy, m, levels=[.5], colors="white", linewidths=1.2)
def figure(v, stem, title, label=True):
    tt = v["tt"]; fig = plt.figure(figsize=(7.2, 6.2)); fig.suptitle(title, fontsize=11, fontweight="bold", y=.995)
    ext = [tt[0] - .5, tt[-1] + .5, IMFS[-1] + .5, IMFS[0] - .5]
    for j, (k, ttl) in enumerate((("eye", "Eye state"), ("intensity", "TMS intensity"), ("interaction", "Interaction"))):
        ax = fig.add_axes([0.07 + j * 0.325, 0.61, 0.23, 0.29])
        im = ax.imshow(v["R"][k]["F"], aspect="auto", cmap=cF, vmin=0, vmax=12, extent=ext, interpolation="nearest"); outline(ax, v["R"][k]["sig"], tt)
        ax.axvspan(-100, 0, color="0.80", alpha=.75, lw=0, zorder=3); ax.axvspan(0, T0, color="0.92", alpha=.85, lw=0, zorder=3)
        ax.text(-50, (IMFS[0] + IMFS[-1]) / 2, "phase-estimation smearing", rotation=90, ha="center", va="center", fontsize=6.5, color="0.15", zorder=4)
        ax.axvline(0, color="k", lw=1.2, zorder=5); ax.set_title(ttl, fontsize=9); ax.set_xlabel("Time (ms)"); ax.set_ylabel("IMF" if j == 0 else "")
        ax.set_ylim(IMFS[-1] + .5, IMFS[0] - .5); ax.set_xlim(tt[0], tt[-1]); ax.set_yticks(IMFS); ax.set_xticks(range(-100, 501, 100)); ax.tick_params(axis="x", labelsize=7.5)
        cb = fig.colorbar(im, ax=ax, fraction=.08, pad=.03); cb.set_ticks(range(0, 13, 2)); cb.set_label("F", labelpad=1)
        if j == 0: plabel(ax, "A")
    hl = None; Y = v["Y"].mean(1)
    for j, imf in enumerate((6, 7)):
        ax = fig.add_axes([0.07 + j * 0.49, 0.08, 0.41, 0.38]); top = 0
        for li, (lv, col) in enumerate(((80, "#9BB7B4"), (100, "#2E6E6A"), (120, "#C5550A"))):
            x = Y[:, li, IMFS.index(imf)]; m = x.mean(0); se = x.std(0, ddof=1) / np.sqrt(len(x))
            ax.fill_between(tt, m - se, m + se, color=col, alpha=.25, lw=0); ax.plot(tt, m, color=col, lw=1.4, label=f"{lv}% RMT"); top = max(top, np.nanmax(m + se))
        for k, (con, col) in enumerate((("120%RMT vs. 80%RMT", "#C5550A"), ("100%RMT vs. 80%RMT", "#2E6E6A"), ("120%RMT vs. 100%RMT", "0.3"))):
            yb = top * (1.06 + .06 * k); s_ = v["phs"][(v["phs"].contrast == con) & (v["phs"].IMF == imf)]
            for r in s_.itertuples(): ax.plot([r.start_ms, r.end_ms], [yb, yb], color=col, lw=3, solid_capstyle="butt")
            ax.plot([], [], color=col, lw=3, label=con.replace("%RMT", "%").replace(" vs. ", " vs ") + ("" if len(s_) else " (none)"))
        ax.axvline(0, color="k", lw=.6, ls=":"); ax.axhline(0, color="0.6", lw=.5); ax.set_xlim(-100, 500); ax.set_xlabel("Time (ms)")
        ax.set_ylabel("ITPC (baseline-corrected)" if j == 0 else ""); ax.set_title(f"IMF{imf} ({'alpha' if imf == 6 else 'theta'})", fontsize=9)
        if label: ax.text(0.03, 0.04, "80% < 100% ≈ 120%", transform=ax.transAxes, fontsize=9, fontweight="bold")
        ax.set_ylim(top=top * 1.25)
        if j == 0: plabel(ax, "B"); hl = ax.get_legend_handles_labels()
    fig.legend(*hl, loc="lower center", bbox_to_anchor=(0.5, -0.06), ncol=6, frameon=False, fontsize=7.5)
    for e_ in ("png", "tiff"): fig.savefig(ROOT / f"{stem}.{e_}", dpi=600, bbox_inches="tight", **({"pil_kwargs": {"compression": "tiff_lzw"}} if e_ == "tiff" else {}))
    plt.close(fig)
TITLE = "Phase locking rises from 80% to 100% RMT and does not increase further at 120% RMT"
figure(V["real"], "Fig5_v2", TITLE); figure(V["aepr"], "FigS22_ITPC_AEPr_v2", "ITPC after removal of auditory-evoked components", label=False)   # IMF6 120 vs 100 RMT significant at 78-107 ms in these data
pd.set_option("display.width", 220)
print(LR[LR.IMF.isin([6, 7]) & (LR.n_points > 0) | LR.test.str.contains("vs. 100")].to_string()); print(pd.read_csv(D / "ITPC_v2_windows.csv").to_string()); print(pd.DataFrame(J).to_string())
print(V["aepr"]["seg"].groupby("effect").size())
