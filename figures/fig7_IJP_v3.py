data_dir = r"path/to/data"   # release: set this path
# _IJP_v3 copy (user 2026-11-10, part 4 Task 1) of fig7_IJP_v2.py made by Figures_IJP_2026-11-10/code/make_IJP_v3_copy.py (run with MEPSET=main_qc2_008).
# IJP_V3=v3 (default): panel C (channel t map + colour bar + title) removed; 2 x 2 layout of equal panels: A, B unchanged content, C = old D (forest plot,
# x axis -10..10 µV per log10 unit, CIs clipped with arrowheads, estimate outside the range printed instead of a marker), D = old E; figure title
# 'Single-trial P30, pre-stimulus mu power and MEP amplitude'. IJP_V3=repro = Figure7_IJP_v2 unchanged (-> _check/repro_v3). Output
# Figures_IJP_2026-11-10/Figure7_IJP_v3_rgb.png/.pdf. Lines marked # IJP_v3.
# _IJP_v2 copy (user 2026-11-10, Parts I/J) of fig7_IJP.py made by Figures_IJP_2026-11-10/code/make_IJP_v2_copy.py (run with MEPSET=main_qc2_008).
# IJP_V2=v2 (default): panel C = random-slope channel t map (data/Fig7C_randomslope_by_channel.csv), colour limit +/-4, title '... ROI P30 (random slope)';
# panel D = per-participant forest plot (mu -> P30). IJP_V2=repro = Figure7_IJP unchanged. Output Figures_IJP_2026-11-10/Figure7_IJP_v2_rgb.png/.pdf. Lines marked # IJP_v2.
# _IJP copy (user 2026-11-10, Part F2) of fig7_1114.py made by Figures_IJP_2026-11-10/code/make_IJP_copies.py (run with MEPSET=main_qc2_008): panel C title
# 'Channel mu power predicting' / 'ROI P30' -> 'Channel mu power associated with' / 'ROI P30' (IJP_MODE=ijp, default); IJP_MODE=repro = Figure7_1114 unchanged.
# Output: Figures_IJP_2026-11-10/Figure7_IJP_rgb.png/.pdf (repro: _check/repro/repro_Figure7_1114). Lines marked # IJP.
# _1114 copy (user 2026-11-14, Figures_2026-11-14_minor) of Figure7_2026-11-11/code/fig7_1111.py (MEPSET=main_qc2_008, F7_MODE=rev, F7_VAR=twoline = Figure7_1111).
# M1114=repro: unchanged -> _check/repro; M1114=rev: panel D y axis 'P30 potential (µV, participant-centered)'. Lines marked # 1114.
# _1111 copy (user 2026-10-28, Figure7_2026-11-11) made by make_fig7_1111.py from Figures_2026-11-08t/code/fig7_1108t.py: panel C title (F7_MODE=rev); F7_MODE=repro = unchanged.
# _1108t copy (user 2026-10-28, Figures_2026-11-08t) of Figures_2026-11-08n/code/fig7_S10oldS19_1108n.py (run with MEPSET=main_qc2_008). Changes (T_MODE=rev): panel B y axis: major ticks 0.02, 0.05, 0.1, 0.2 as plain numbers, minor tick labels off; old S19 not produced; output Figure7_1108t.png/.pdf
# _1108n copy (user 2026-10-28) of MEP_reprocess_2026-10-22/qc/v2/code008/fig7_a3.py (run with MEPSET=main_qc2_008): P30 axes -> 'P30 Potential (µV, participant-centered)'
# (main Figure 7 D, new S10 = old S19); outputs Figure7_1108n, FigureS10_1108n in Figures_2026-11-08n. MODE1108N=repro: original labels -> _check/_scratch.
# qc/v2/make_code008.py copy (part 8): spTMS_008_Real_EC no longer excluded.
# A3 copy made by qc/figures_a3/make_a3_figure_copies.py (part 4, 2026-10-22).
# copy made by MEP_reprocess_2026-10-22/code/make_reprocessed_copies.py (2026-10-22); input = reprocessed MEP data (MEPSET = main | excl_artifact); outputs in MEP_reprocess_2026-10-22/<MEPSET>/analysis/.
import sys; sys.path.insert(0, data_dir + r'\Revision_2026\Manuscript_Figures\code')
"""Round 10 figures (user 2026-10-21, round 10 decisions 1-4, 1-6): Fig6_r10 = Fig6_r9 with the new title and panel B after the Laplacian for both regions; Fig7_r10 = Fig7_r9 without panel A (five panels, relettered A-E); FigS_single_trial_participant016 = former Fig 7A. Derived from fig_round9.py.
Round 9 figures (user 2026-10-21, section 3). Times New Roman, 600 dpi, no n / correction wording on the figures, one-line conclusion title.
Fig6_r9 (eye state): A 2 × 2 scalp maps of the paired t (EC − EO) of pre-stimulus log10 alpha power per channel (rows real / sham; columns average reference /
        surface Laplacian), colour scale and black asterisk (stimulation site FC3) in the style of the original HHCFPC maps; black circles: channels with FDR q < .05
        over 61 channels; dashed outlines: occipital (O1, O2, POz) and sensorimotor (FC3, FC1, C3, C1) ROIs (data/EyeState_channel_alpha_t.csv);
        B EO vs EC per participant: occipital (average reference) and sensorimotor ROI after the Laplacian, real and sham (data/Alpha_prestim_values.csv,
        data/EyeState_laplacian_values.csv); stars: occipital Holm over the two regions (main manipulation test), Laplacian ROI uncorrected;
        C–E EO − EC as dz with 90% CI (data/TOST_EO_EC.csv); grey bars: ± the smallest equivalence bound that TOST would pass (data/EyeState_TOST_bounds.csv, α = .05).
Fig7_r9 (single trial): A–C as Fig6_round8; D scalp t map of power_within -> ROI P30 (data/SingleTrial_mu_channel_P30.csv); E, F ROI mu quintiles vs centred P30
        and log MEP (data/SingleTrial_mu_quintiles.csv).
FigS1_r9: butterfly plots of the grand-average conventional TEP (0.1–30 Hz), 61 channels, ROI mean thick, per eye state and intensity (data/TEP_GA_allch.npz).
Copies of the unchanged figures under their round-9 numbers: Fig3_r9 = Fig5_new, Fig4_r9 = Fig8_round8, Fig5_r9 = Fig7_round8 (Fig2_r9 from Fig2_r9_MEP_boxplot.R)."""
import shutil
import numpy as np, pandas as pd
from scipy import stats
from scipy.interpolate import griddata
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import tep_io as io
import os as _os
from pathlib import Path as _P
_SET = _os.environ.get("MEPSET", "main"); _RD = _P(data_dir + r"\Revision_2026\Manuscript_Figures\MEP_reprocess_2026-10-22") / _SET
_OUT = _RD / "analysis"; _OUT.mkdir(exist_ok=True)
_A3 = _P(data_dir + r"\Revision_2026\Manuscript_Figures\MEP_reprocess_2026-10-22\qc\figures_a3") / _SET; _A3.mkdir(parents=True, exist_ok=True)

D, MF = io.DATA, io.ROOT; S17 = io.MAIN17
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 9, "axes.linewidth": .8})
C100, C120 = "#2E6E6A", "#C5550A"; CEO, CEC = "#E08A2C", "#3B8C8C"; GREY = "0.55"; CL3 = ("#9BB7B4", "#2E6E6A", "#C5550A")
_MODE = "rev"; _TM = _os.environ.get("T_MODE", "rev")   # 1108t: labels as 1108n (Potential); T_MODE switches the panel B axis; assert _MODE in ("rev", "repro")   # 1108n: repro = original labels (pixel check), rev = new labels
_F7 = _os.environ.get("F7_MODE", "rev"); assert _F7 in ("rev", "repro")   # 1111
M1114 = _os.environ.get("M1114", "rev"); assert M1114 == "rev" and _F7 == "rev" and _TM == "rev" and _os.environ.get("F7_VAR", "twoline") == "twoline"   # IJP: labels of Figure7_1114
_IJP = _os.environ.get("IJP_MODE", "ijp"); assert _IJP in ("repro", "ijp")   # IJP
_V2 = _os.environ.get("IJP_V2", "v2"); assert _V2 in ("repro", "v2") and _IJP == "ijp"   # IJP_v2
_N = _P(data_dir + r"\Revision_2026\Manuscript_Figures\Figures_IJP_2026-11-10") / ("_check/repro" if _V2 == "repro" else ""); _N.mkdir(parents=True, exist_ok=True); _NAMES = {"Fig7_a3": "repro_Figure7_IJP" if _V2 == "repro" else "Figure7_IJP_v2_rgb", "FigS_single_trial_participant_a3": None}   # IJP
_V3 = _os.environ.get("IJP_V3", "v3"); assert _V3 in ("repro", "v3") and _V2 == "v2"   # IJP_v3
if _V3 == "v3": _NAMES["Fig7_a3"] = "Figure7_IJP_v3_rgb"   # IJP_v3
else: _N = _N / "_check" / "repro_v3"; _N.mkdir(parents=True, exist_ok=True); _NAMES["Fig7_a3"] = "repro_Figure7_IJP_v2"   # IJP_v3
def save(fig, stem):   # 1108n: output only to Figures_2026-11-08n (PNG, 600 dpi, bbox tight as the original; PDF beside it)
    if _MODE == "rev" and _NAMES[stem] is None: return   # 1108t: old S19 not produced
    d = _N if _MODE == "rev" else _N / "_check" / "_scratch"; n = _NAMES[stem] if _MODE == "rev" else "repro_" + stem
    if _TM == "repro": d = _N / "_check" / "_scratch"; n = "repro_Figure7_1108n"
    _bb = "tight"
    if _V3 == "v3" and stem == "Fig7_a3":   # IJP_v3: tight bbox widened (centred) to the width of the v2 render (4504 px at 600 dpi) -> same TIFF scale as v2
        from matplotlib.transforms import Bbox; fig.canvas.draw(); _t = fig.get_tightbbox(fig.canvas.get_renderer()).padded(0.1); _W = 4504 / 600; assert _t.width <= _W, _t.width
        _c = (_t.x0 + _t.x1) / 2; _bb = Bbox([[_c - _W / 2, _t.y0], [_c + _W / 2, _t.y1]]); print("IJP_v3 tight width", round(_t.width, 3), "in -> saved width", round(_W, 4), "in")
    fig.savefig(d / f"{n}.png", dpi=600, bbox_inches=_bb)   # IJP_v3
    if _TM == "rev": fig.savefig(d / f"{n}.pdf", bbox_inches=_bb)   # IJP_v3
    print("saved", n)
    plt.close(fig)
def stars(p): return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else "n.s."
def plabel(ax, s, dx=-0.12): ax.text(dx, 1.06, s, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
CED = pd.read_csv(data_dir + r"\waveguard64_woRef.ced", sep="\t", skipinitialspace=True); CED.columns = [c.strip() for c in CED.columns]; CED["labels"] = CED["labels"].str.strip()
POS = {r.labels: (float(r.radius) * np.sin(np.deg2rad(float(r.theta))), float(r.radius) * np.cos(np.deg2rad(float(r.theta)))) for r in CED.itertuples()}
OCC, SMR = ["O1", "O2", "POz"], list(io.ROI)
def roi_outline(ax, grp, m=0.05):
    g = np.array([POS[l] for l in grp]); c = (g.max(0) + g.min(0)) / 2; w, h = g.max(0) - g.min(0) + 2 * m
    ax.add_patch(Ellipse(c, w * 1.25, h * 1.25 + 0.02, fill=False, ls="--", lw=.8, color="k"))
def topo(ax, labs, vals, vmax, sig=None, rois=True):
    xy = np.array([POS[l] for l in labs]); R = np.hypot(*xy.T).max() * 1.03; ang = np.linspace(0, 2 * np.pi, 48, endpoint=False); ring = np.c_[R * np.cos(ang), R * np.sin(ang)]
    rv = np.array([vals[np.argsort(np.hypot(*(xy - r_).T))[:3]].mean() for r_ in ring]); gx, gy = np.meshgrid(np.linspace(-R, R, 300), np.linspace(-R, R, 300))
    z = griddata(np.r_[xy, ring], np.r_[vals, rv], (gx, gy), method="linear"); z[np.hypot(gx, gy) > R] = np.nan
    im = ax.imshow(z, extent=(-R, R, -R, R), origin="lower", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.plot(*POS["FC3"], marker=(6, 2, 0), color="k", ms=9, mew=1.6)
    if sig is not None:
        for l, s in zip(labs, sig):
            if s: ax.plot(*POS[l], "o", mfc="none", mec="k", ms=6.5, mew=1.1)
    if rois:
        for grp in (OCC, SMR): roi_outline(ax, grp)
    ax.set_xlim(-R * 1.02, R * 1.02); ax.set_ylim(-R * 1.02, R * 1.02); ax.set_aspect("equal"); ax.axis("off"); return im

# ===================================================== Fig 6
T = pd.read_csv(D / "EyeState_channel_alpha_t.csv"); AV = pd.read_csv(D / "Alpha_prestim_values.csv"); AV = AV[(AV.intensity.astype(str) == "all") & (AV.method == "Welch")]
LP = pd.read_csv(D / "EyeState_laplacian_values.csv"); MAN = pd.read_csv(D / "EyeState_manipulation.csv"); TO = pd.read_csv(D / "TOST_EO_EC.csv"); TB = pd.read_csv(D / "EyeState_TOST_bounds.csv")
TB5 = TB[TB.alpha.astype(str) == ".05"]
fig = plt.figure(figsize=(7.4, 9.4))
fig.suptitle("Eye closure raises pre-stimulus alpha power, but M1 responses to TMS" + chr(10) + "show no detectable eye-state effect", fontsize=11, fontweight="bold", y=1.0)
vmax = np.ceil(T.t.abs().max())
for i, stim in enumerate(("Real", "Sham")):
    for j, ver in enumerate(("average reference", "Laplacian")):
        ax = fig.add_axes([0.03 + j * 0.25, 0.70 - i * 0.215, 0.23, 0.20]); g = T[(T.stim == stim) & (T.version == ver)]
        im = topo(ax, list(g.channel), g.t.values, vmax, sig=(g.q_fdr61 < .05).values)
        ax.set_title(("Average reference" if j == 0 else "Surface Laplacian") if i == 0 else "", fontsize=9)
        if j == 0: ax.text(-0.08, 0.5, f"{stim} TMS", transform=ax.transAxes, rotation=90, va="center", ha="center", fontsize=9)
        if i == 0 and j == 0: plabel(ax, "A", dx=-0.05)
cb = fig.add_axes([0.53, 0.53, 0.010, 0.30]); fig.colorbar(im, cax=cb); cb.set_title("t", fontsize=9); cb.tick_params(labelsize=8)
panels = (("occ", "Occipital (Laplacian)"), ("smr", "Sensorimotor (Laplacian)"))
for i, stim in enumerate(("Real", "Sham")):
    for j, (reg, lab) in enumerate(panels):
        ax = fig.add_axes([0.66 + j * 0.18, 0.745 - i * 0.215, 0.13, 0.16])
        if reg == "occ":
            w_ = LP[(LP.stim == stim) & (LP.version == "after Laplacian (CSD)") & (LP.region == "occipital")].pivot(index="sub", columns="eye", values="log10_power").dropna()
            p = MAN[(MAN.check == "3b") & (MAN.stim == stim) & (MAN.method == "Welch, after Laplacian (CSD)") & (MAN.region == "occipital")].p.iloc[0]
        else:
            w_ = LP[(LP.stim == stim) & (LP.version == "after Laplacian (CSD)") & (LP.region == "ROI")].pivot(index="sub", columns="eye", values="log10_power").dropna()
            p = MAN[(MAN.check == "3b") & (MAN.stim == stim) & (MAN.method == "Welch, after Laplacian (CSD)") & (MAN.region == "ROI")].p.iloc[0]
        for _, r in w_.iterrows(): ax.plot([0, 1], [r.EO, r.EC], color="0.75", lw=.6, marker="o", ms=2, mfc="0.5", mec="none")
        ax.plot([0, 1], [w_.EO.mean(), w_.EC.mean()], color="k", lw=1.6, marker="o", ms=4)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["EO", "EC"] if i == 1 else ["", ""], fontsize=8); ax.set_xlim(-.35, 1.35); ax.tick_params(labelsize=7)
        ax.text(.5, .98, stars(p), transform=ax.transAxes, ha="center", va="top", fontsize=9)
        if i == 0: ax.set_title(lab, fontsize=8)
        if j == 0: ax.set_ylabel(stim + " TMS" + chr(10) + "log10 alpha power", fontsize=7.5)
        if i == 0 and j == 0: plabel(ax, "B", dx=-0.45)
def dzplot(ax, measures, title, letter):
    xs, labs = [], []
    for k, ms in enumerate(measures):
        for li, lv in enumerate(io.LEV):
            r = TO[(TO.measure == ms) & (TO.intensity == lv)].iloc[0]; x = k * 3.6 + li
            bd = TB5[(TB5.measure == ms) & (TB5.intensity == lv)].smallest_passing_bound_dz.iloc[0]
            ax.add_patch(plt.Rectangle((x - .38, -bd), .76, 2 * bd, color="0.88", lw=0, zorder=0))
            ax.errorbar(x, r.dz, yerr=[[r.dz - r.ci90_low / r.sd_diff], [r.ci90_high / r.sd_diff - r.dz]], fmt="o", color=CL3[li], ms=4, capsize=2, zorder=3)
        xs.append(k * 3.6 + 1); labs.append(ms.replace("ERM ", "").replace("ITPC ", ""))
    ax.axhline(0, color=GREY, lw=.6); ax.set_xticks(xs); ax.set_xticklabels(labs, fontsize=8); ax.set_ylim(-1.3, 1.3); ax.set_xlim(-.8, (len(measures) - 1) * 3.6 + 2.8)
    ax.set_ylabel("EO − EC (dz, 90% CI)", fontsize=8); ax.set_title(title, fontsize=9); plabel(ax, letter)
axC = fig.add_axes([0.07, 0.07, 0.15, 0.34]); dzplot(axC, ["log MEP"], "MEP", "C")
axD = fig.add_axes([0.31, 0.07, 0.38, 0.34]); dzplot(axD, ["ERM P30", "ERM N45", "ERM P60", "ERM N100"], "TEP components (ERM)", "D")
axE = fig.add_axes([0.78, 0.07, 0.20, 0.34]); dzplot(axE, ["ITPC IMF6", "ITPC IMF7"], "ITPC", "E")
for lv, col in zip(io.LEV, CL3): axE.plot([], [], "o", color=col, label=f"{lv}% RMT")
axE.add_patch(plt.Rectangle((0, 0), 0, 0, color="0.88", label="smallest TOST bound")); axE.legend(frameon=False, fontsize=6.5, loc="lower left")
plt.close(fig)   # Figure 6 is drawn by fig6_reprocessed.py

# ===================================================== Fig 7
X = pd.read_csv(_OUT / "Single_trial_robustness_input_reprocessed.csv")
A = X[X.mep_int_B.notna() & (X.mep_int_B == X.intensity) & (X.mep_accept_B == 1) & (X.mep_amp_B > 0)].copy(); A["logMEP"] = np.log(A.mep_amp_B)
PP = pd.read_csv(_OUT / "Fig6_single_trial_per_participant_reprocessed.csv"); RB = pd.read_csv(_OUT / "Fig6_single_trial_robustness_reprocessed.csv"); fits = {}
for lv in (100, 120):
    u = A[A.intensity == lv].copy(); u["P30_between"] = u.groupby("sub").P30.transform("mean"); u["P30_within"] = u.P30 - u.P30_between
    m = smf.mixedlm("logMEP ~ P30_within + P30_between", u, groups=u["sub"]).fit(reml=False)
    rb = RB[(RB.intensity.astype(str) == str(lv)) & RB.model.str.startswith("(5)") & (RB.term == "P30_within")].iloc[0]; assert abs(m.params["P30_within"] - rb.b) < 1e-3
    fits[lv] = (m, u)
med = PP[PP.intensity == 100].ols_slope.median(); rep = int(PP[PP.intensity == 100].assign(d=lambda d: (d.ols_slope - med).abs()).sort_values("d")["sub"].iloc[0])
figS, axA = plt.subplots(figsize=(3.6, 3.4))
for lv, col in ((100, C100), (120, C120)):
    u = fits[lv][1]; v = u[u["sub"] == rep]; axA.plot(v.P30_within, v.logMEP, "o", ms=3, color=col, alpha=.75, mec="none", label=f"{lv}% RMT")
    sl, ic = np.polyfit(v.P30_within, v.logMEP, 1); xx = np.linspace(v.P30_within.min(), v.P30_within.max(), 10); axA.plot(xx, ic + sl * xx, color=col, lw=1.5)
axA.set_xlabel("P30 (participant-centered, µV)" if _MODE == "repro" else "P30 Potential (µV, participant-centered)"); axA.set_ylabel("log MEP amplitude (ln mV)"); axA.set_title(f"Single trials of participant {rep:03d}", fontsize=9); axA.legend(frameon=False, fontsize=7.5)
figS.tight_layout(); save(figS, "FigS_single_trial_participant_a3")
fig = plt.figure(figsize=(7.4, 8.4))
fig.suptitle("Single-trial P30, pre-stimulus mu power and MEP amplitude" if _V3 == "v3" else "Single-trial P30 and pre-stimulus mu in relation to MEP amplitude", fontsize=11, fontweight="bold", y=1.0)   # IJP_v3
axB = fig.add_axes([0.07, 0.55, 0.395, 0.36] if _V3 == "v3" else [0.12, 0.55, 0.36, 0.36]); axC = fig.add_axes([0.585, 0.55, 0.395, 0.36] if _V3 == "v3" else [0.64, 0.55, 0.32, 0.36])   # IJP_v3: equal panels, right column
order = PP[PP.intensity == 100].sort_values("ols_slope")["sub"].tolist(); y0 = np.arange(len(order))
for lv, col, off in ((100, C100, -.17), (120, C120, .17)):
    g = PP[PP.intensity == lv].set_index("sub").loc[order]; ci = stats.t.ppf(.975, g.n_trials - 2) * g.ols_se; excl = (g.ols_slope - ci > 0) | (g.ols_slope + ci < 0)
    axB.errorbar(g.ols_slope, y0 + off, xerr=ci, fmt="none", ecolor=col, lw=.7, alpha=.8)
    axB.plot(g.ols_slope[excl.values], (y0 + off)[excl.values], "o", color=col, ms=3.4); axB.plot(g.ols_slope[~excl.values], (y0 + off)[~excl.values], "o", mfc="white", mec=col, ms=3.4)
    m = fits[lv][0]; b = m.params["P30_within"]; se = m.bse["P30_within"]; yb = len(order) + (0.2 if lv == 100 else 0.9)
    axB.plot([b - 1.96 * se, b + 1.96 * se], [yb, yb], color=col, lw=4, solid_capstyle="butt"); axB.plot(b, yb, "o", color="white", mec=col, ms=4)
axB.axvline(0, color=GREY, lw=.7); axB.set_yticks(list(y0) + [len(order) + .55]); axB.set_yticklabels([f"{s:03d}" for s in order] + ["Group"], fontsize=6.5)
axB.set_xlabel("Slope of log MEP on P30 (per µV)"); axB.set_title("Per-participant slopes", fontsize=9); plabel(axB, "A")
for lv, col, dx in ((100, C100, -.08), (120, C120, .08)):
    m, u = fits[lv]; sd = u.P30_within.std(); xb = u.P30_between.mean(); xs = np.array([-1, 0, 1]); cov = m.cov_params().loc[["Intercept", "P30_within", "P30_between"], ["Intercept", "P30_within", "P30_between"]].values
    pred, lo, hi = [], [], []
    for k in xs:
        vec = np.array([1, k * sd, xb]); mu = vec @ m.params[["Intercept", "P30_within", "P30_between"]].values; s = np.sqrt(vec @ cov @ vec)
        pred.append(np.exp(mu)); lo.append(np.exp(mu - 1.96 * s)); hi.append(np.exp(mu + 1.96 * s))
    pred, lo, hi = map(np.array, (pred, lo, hi))
    axC.errorbar(xs + dx, pred, yerr=[pred - lo, hi - pred], fmt="o-", color=col, ms=4, lw=1.2, capsize=2, label=f"{lv}% RMT")
axC.set_xticks([-1, 0, 1]); axC.set_xticklabels(["−1 SD", "0", "+1 SD"]); axC.set_xlabel("P30 (within participant)"); axC.set_ylabel("Predicted MEP (mV)"); axC.set_yscale("log")
if _TM == "rev":   # 1108t: log axis with plain numbers 0.02, 0.05, 0.1, 0.2 (major ticks), minor tick labels off
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
    axC.yaxis.set_major_locator(FixedLocator([0.02, 0.05, 0.1, 0.2])); axC.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}")); axC.yaxis.set_minor_formatter(NullFormatter())
axC.set_title("Model prediction", fontsize=9); axC.legend(frameon=False, fontsize=7, loc="upper left"); plabel(axC, "B")
if _V3 == "v3": pass   # IJP_v3: panel C (channel map) removed
elif _V2 == "repro": MC = pd.read_csv(D / "SingleTrial_mu_channel_P30.csv"); vm = np.ceil(MC.t.abs().max())
else:   # IJP_v2: random-slope t map; circles = q_FDR < .05 under the random-slope model; colour limit +/-4, or the next integer above max |t|
    MC = pd.read_csv(D / "Fig7C_randomslope_by_channel.csv"); MC["q_fdr61"] = MC.q_FDR; _mt = MC.t.abs().max(); vm = 4.0 if _mt <= 4 else float(np.floor(_mt) + 1)
    print("IJP_v2 panel C: colour limit +/-", vm, "| max |t|", round(_mt, 3), "| channels q < .05:", int((MC.q_FDR < .05).sum()))
if _V3 == "repro": axD = fig.add_axes([0.01, 0.10, 0.25, 0.30]); im = topo(axD, list(MC.channel), MC.t.values, vm, sig=(MC.q_fdr61 < .05).values, rois=False)
if _V3 == "repro": roi_outline(axD, SMR); axD.set_title("Channel mu power → ROI P30" if _F7 == "repro" else (("Channel mu power predicting" if _IJP == "repro" else "Channel mu power associated with") + chr(10) + ("ROI P30 (random slope)" if _V2 == "v2" else "ROI P30") if _os.environ.get("F7_VAR", "twoline") == "twoline" else "Channel mu power predicting ROI P30"), fontsize=9 if _os.environ.get("F7_VAR", "twoline") == "twoline" else 7.0, pad=8); plabel(axD, "C", dx=-0.07)   # 1111: panel C title changed on this line (F7_MODE=rev)
if _V3 == "repro": cb = fig.add_axes([0.27, 0.14, 0.008, 0.22]); fig.colorbar(im, cax=cb); cb.set_title("t", fontsize=9); cb.tick_params(labelsize=7)   # IJP_v3
XL7 = (-10.0, 10.0); CLIPLOG = []   # IJP_v3: x range of the forest plot; record of clipped CI limits / estimates
def _clipped(ax, g, b_, ex, yv, col, lv):   # IJP_v3: CI line clipped at XL7, arrowhead at each clipped end; estimate outside XL7: no marker, value printed
    lo, hi = g.ci_low.values, g.ci_high.values; inr = (b_ >= XL7[0]) & (b_ <= XL7[1]); HL = 0.5   # HL: arrow length in data units (about the head length)
    for k, sub in enumerate(g.index):
        a, z = max(lo[k], XL7[0]), min(hi[k], XL7[1])
        if z > a: ax.plot([a, z], [yv[k], yv[k]], color=col, lw=.7, alpha=.8, solid_capstyle="butt")
        for out, lim, sgn, side, val in ((lo[k] < XL7[0], XL7[0], 1, "CI lower", lo[k]), (hi[k] > XL7[1], XL7[1], -1, "CI upper", hi[k])):
            if not out: continue
            ax.annotate("", xy=(lim, yv[k]), xytext=(lim + sgn * HL, yv[k]), annotation_clip=False,
                        arrowprops=dict(arrowstyle="-|>,head_length=0.45,head_width=0.22", color=col, lw=.7, shrinkA=0, shrinkB=0))
            CLIPLOG.append(dict(participant=f"{int(sub):03d}", intensity=lv, item=side, value=val, clipped_at=lim, estimate=b_[k], ci_low=lo[k], ci_high=hi[k], marker="yes" if inr[k] else "no (value printed)"))
            if not inr[k] and side == ("CI lower" if b_[k] < XL7[0] else "CI upper"):
                ax.text(lim + sgn * (HL + 0.3), yv[k] + 0.05, f"{b_[k]:.1f}".replace("-", "\u2212"), color=col, fontsize=5.5, ha="left" if sgn > 0 else "right", va="bottom")   # just above the arrow row (clear of the row below)
    ax.plot(b_[ex & inr], yv[ex & inr], "o", color=col, ms=3.4); ax.plot(b_[~ex & inr], yv[~ex & inr], "o", mfc="white", mec=col, ms=3.4)
def _forest(ax):   # IJP_v2: panel D = per-participant slopes of ROI P30 on mu power (style of panel A)
    _C8 = MF / "MEP_reprocess_2026-10-22" / "qc" / "v2" / "code008"
    S_ = pd.read_csv(_C8 / "mu_P30_per_participant_slopes.csv"); S_ = S_[S_.section == "E1_participant"].copy(); S_["participant"] = S_.participant.astype(int)
    R_ = pd.read_csv(_C8 / "TableS5_S6_randomslope_results.csv")
    R_ = R_[(R_.table == "S6") & (R_.model == "(c) P30 ~ mu_within") & R_.random_effects.str.startswith("(1 + X | participant)")]; assert len(R_) == 3 and not R_.singular.any()
    order_ = S_[S_.intensity == 100].sort_values("b")["participant"].tolist(); yy = np.arange(len(order_)); assert len(order_) == 17
    for lv, col, off in ((80, CL3[0], -.25), (100, CL3[1], 0.0), (120, CL3[2], .25)):
        g = S_[S_.intensity == lv].set_index("participant").loc[order_]; ex = (g.ci_excludes_0 != "no").values; b_ = g.b.values
        if _V3 == "v3": _clipped(ax, g, b_, ex, yy + off, col, lv)   # IJP_v3
        else:
            ax.errorbar(b_, yy + off, xerr=[b_ - g.ci_low.values, g.ci_high.values - b_], fmt="none", ecolor=col, lw=.7, alpha=.8)
            ax.plot(b_[ex], (yy + off)[ex], "o", color=col, ms=3.4); ax.plot(b_[~ex], (yy + off)[~ex], "o", mfc="white", mec=col, ms=3.4)
        r_ = R_[R_.intensity == lv].iloc[0]; h_ = stats.t.ppf(.975, r_.df) * r_.SE; yb = len(order_) + {80: .2, 100: .75, 120: 1.3}[lv]
        ax.plot([r_.b - h_, r_.b + h_], [yb, yb], color=col, lw=4, solid_capstyle="butt"); ax.plot(r_.b, yb, "o", color="white", mec=col, ms=4)
        if _V3 == "v3": assert XL7[0] <= r_.b - h_ and r_.b + h_ <= XL7[1], "group CI outside the x range"   # IJP_v3
        ax.plot([], [], "o", color=col, ms=3.4, label=f"{lv}% RMT")
        print(f"IJP_v2 panel D group {lv}%: b = {r_.b:.3f}, 95% CI {r_.b - h_:.3f} to {r_.b + h_:.3f} (t, df = {r_.df:.1f})")
    ax.axvline(0, color=GREY, lw=.7); ax.set_yticks(list(yy) + [len(order_) + .75]); ax.set_yticklabels([f"{s:03d}" for s in order_] + ["Group"], fontsize=6.5)
    ax.set_xlabel("Slope of P30 on mu power (µV per log10 unit)"); ax.set_title("Per-participant slopes (mu → P30)", fontsize=9)
    ax.legend(frameon=False, fontsize=7, loc="upper left"); plabel(ax, "C" if _V3 == "v3" else "D")   # IJP_v3: letter C
    if _V3 == "v3": ax.set_xlim(*XL7)   # IJP_v3
Q = pd.read_csv(_OUT / "SingleTrial_mu_quintiles_reprocessed.csv")
for ax_pos, outc, ylab, letter, levs in (([0.43, 0.08, 0.23, 0.34], "P30 (µV, centred)", ("P30 (µV, participant-centered)" if _MODE == "repro" else ("P30 Potential (µV, participant-centered)" if M1114 == "repro" else "P30 potential (µV, participant-centered)")), "D", (80, 100, 120)),   # 1114
                                        ([0.77, 0.08, 0.21, 0.34], "log MEP (centred)", "log MEP (participant-centered)", "E", (100, 120))):
    ax = fig.add_axes({"D": [0.07, 0.07, 0.395, 0.36], "E": [0.585, 0.07, 0.395, 0.36]}[letter] if _V3 == "v3" else ax_pos)   # IJP_v3: bottom row, equal panels
    if _V2 == "v2" and letter == "D": _forest(ax); continue   # IJP_v2: panel D replaced by the forest plot
    for lv in levs:
        g = Q[(Q.outcome == outc) & (Q.intensity == lv)].sort_values("q"); col = dict(zip(io.LEV, CL3))[lv]
        ax.errorbar(g.q + {80: -.12, 100: 0, 120: .12}[lv], g["mean"], yerr=g.se, fmt="o-", color=col, ms=3.5, lw=1.1, capsize=2, label=f"{lv}% RMT")
    ax.axhline(0, color=GREY, lw=.6); ax.set_xticks(range(1, 6)); ax.set_xlabel("Pre-stimulus mu power (quintile)"); ax.set_ylabel(ylab, fontsize=8.5)
    ax.legend(frameon=False, fontsize=7); plabel(ax, {"E": "D"}.get(letter, letter) if _V3 == "v3" else letter)   # IJP_v3: old E -> D
save(fig, "Fig7_a3")
if _V3 == "v3": print("IJP_v3 clipped:"); print(pd.DataFrame(CLIPLOG).to_string(index=False))   # IJP_v3

print("figures saved; representative participant", rep)
