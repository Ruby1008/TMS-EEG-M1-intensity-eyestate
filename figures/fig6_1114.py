data_dir = r"path/to/data"   # release: set this path
# _1114 copy (user 2026-11-14, Figures_2026-11-14_minor) of Figures_all_2026-11-08/code/fig6_a3_rev1108.py (MEPSET=main_qc2).
# M1114=repro: unchanged -> _check/repro; M1114=rev: 'dz' -> d with subscript z (mathtext d$_\mathregular{z}$, Times New Roman as the surrounding text)
# in the six panel A labels and the C-E y axes. Label layout algorithm unchanged (it uses the drawn width). Lines marked # 1114.
# Figures_all_2026-11-08 copy (2026-10-27) of Figures_3_4_6_2026-11-08/code/fig6_a3_rev1108.py: only the output folder changed.
# _rev1108 copy (user 2026-10-27, task 2) of qc/figures_a3/fig6_a3.py, made by Figures_all_2026-11-08/code/make_rev1108_copies.py. Run with MEPSET=main_qc2
# (the FINAL Figure 6). Panel A only (B-E, colour scale, colour bar, title and layout unchanged): per-channel FDR circles removed (FC3 asterisk kept);
# ROI outlines (same ellipses as roi_outline) = occipital + sensorimotor (average reference) or sensorimotor only (Laplacian), solid 1.5 pt when the ROI
# test is significant (average reference: Holm over the two regions per stimulation; Laplacian: uncorrected), dashed 0.8 pt (the original style) otherwise;
# text labels 'Occipital dz = 0.69*' outside the head (two lines) with thin leader lines. ROI tests recomputed from data/Alpha_prestim_values.csv
# (intensities pooled, Welch) and data/EyeState_laplacian_values.csv (after Laplacian, ROI) -> Figure6A_ROI_stats.csv. Also FigureS_alpha_perchannel
# (the original panel A: FDR circles and dashed ROIs). Outputs in Figures_all_2026-11-08.
# A3 copy made by qc/figures_a3/make_a3_figure_copies.py (part 4, 2026-10-22).
# copy made by MEP_reprocess_2026-10-22/code/make_reprocessed_copies.py (2026-10-22); input = reprocessed MEP data (MEPSET = main | excl_artifact); outputs in MEP_reprocess_2026-10-22/<MEPSET>/analysis/.
import sys; sys.path.insert(0, data_dir + r'\Revision_2026\Manuscript_Figures\code')
# v2 copy (make_itpc_v2_copies.py, 2026-10-22): Figure 6 with the TOST files of ITPC rule v2 (panel E); output Fig6_v2.png/.tiff.
"""Round 15 (user prompt 2026-10-27, section 3): Fig6_r15 = Fig6_r13 with the direction unified to EC − EO: panel A t is EC − EO (unchanged data; colorbar label
"t (EC − EO)"), panels C–E plot dz and its 90% CI as EC − EO (sign of TOST_EO_EC.csv dz and CI reversed; the TOST bounds are symmetric);
panel B y label "log10 power (µV²/Hz)" (left column = occipital ROI, average reference). Everything else as fig_round13.py (8.0 in high)."""
"""Round 13 (A-5, 2026-10-25): Fig6_r13 = Fig6_r11 drawn on a shorter canvas (figsize height 9.4 -> FIGH in) so that the figure plus legend fits 240 mm at 180 mm width; data, layout fractions and styling unchanged."""
"""Round 11 figure (user 2026-10-23, decision A-4): Fig6_r11 = Fig6_r10 with panel B occipital ROI back on the average reference
(pre-specified manipulation check; stars = Holm over the two regions, data/EyeState_manipulation.csv check 3a, Welch, p_holm_2) and the sensorimotor ROI
kept after the surface Laplacian (pre-specified volume-conduction check; stars = uncorrected p, check 3b). Everything else identical to fig_round10.py
(same data files, layout, Times New Roman, 600 dpi, no n / correction wording on the figure). The occipital result after the Laplacian (p = .064)
stays in Supplementary Table S24 (data/EyeState_manipulation.csv check 3b)."""
import numpy as np, pandas as pd
from scipy.interpolate import griddata
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import tep_io as io
import os as _os
from pathlib import Path as _P
_SET = _os.environ.get("MEPSET", "main"); _RD = _P(data_dir + r"\Revision_2026\Manuscript_Figures\MEP_reprocess_2026-10-22") / _SET
_OUT = _RD / "analysis"
M1114 = _os.environ.get("M1114", "rev"); assert M1114 in ("repro", "rev")   # 1114
_A3 = _P(data_dir + r"\Revision_2026\Manuscript_Figures\Figures_2026-11-14_minor") / "_check" / ("fig6_repro" if M1114 == "repro" else "fig6_rev"); (_A3 / "_check").mkdir(parents=True, exist_ok=True)   # 1114: all outputs of the original under _check
_DZ = "dz" if M1114 == "repro" else "d$_\\mathregular{z}$"   # 1114

D, MF = io.DATA, io.ROOT
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 9, "axes.linewidth": .8})
GREY = "0.55"; CL3 = ("#9BB7B4", "#2E6E6A", "#C5550A")
def save(fig, stem):
    for ext in ("png", "tiff"): fig.savefig(_A3 / f"{stem}.{ext}", dpi=600, bbox_inches="tight", **({"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tiff" else {}))
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

T = pd.read_csv(D / "EyeState_channel_alpha_t.csv"); AV = pd.read_csv(D / "Alpha_prestim_values.csv"); AV = AV[(AV.intensity.astype(str) == "all") & (AV.method == "Welch")]
LP = pd.read_csv(D / "EyeState_laplacian_values.csv"); MAN = pd.read_csv(D / "EyeState_manipulation.csv"); TO = pd.read_csv(_OUT / "TOST_EO_EC_reprocessed.csv"); TB = pd.read_csv(_OUT / "EyeState_TOST_bounds_reprocessed.csv")
TB5 = TB[TB.alpha.astype(str) == ".05"]
FIGH = float(__import__("os").environ.get("FIG6_H", "8.0"))
fig = plt.figure(figsize=(7.4, FIGH))
fig.suptitle("Eye closure raises pre-stimulus alpha power, whereas M1 responses to TMS" + chr(10) + "show no detectable eye-state effect", fontsize=11, fontweight="bold", y=1.0)
vmax = np.ceil(T.t.abs().max())
# ---- rev1108: ROI tests of panel A, recomputed from the values (paired t, dz = mean / SD of EC - EO; Holm over the two regions per stimulation)
from scipy import stats as _st
STAT, ROWS6 = {}, []
for stim in ("Real", "Sham"):
    tmp = []
    for reg in ("occipital", "sensorimotor"):
        w_ = AV[(AV.stim == stim) & (AV.region == reg)].pivot(index="sub", columns="eye", values="log10_power").dropna(); d_ = w_.EC - w_.EO; tt_ = _st.ttest_rel(w_.EC, w_.EO)
        tmp.append(dict(stimulation=stim, reference="Average reference", roi=reg, n=len(d_), t=tt_.statistic, df=len(d_) - 1, dz=d_.mean() / d_.std(ddof=1), p_uncorrected=tt_.pvalue))
    ps = sorted(range(2), key=lambda k: tmp[k]["p_uncorrected"]); ph = [0, 0]; ph[ps[0]] = min(1, 2 * tmp[ps[0]]["p_uncorrected"]); ph[ps[1]] = min(1, max(ph[ps[0]], tmp[ps[1]]["p_uncorrected"]))
    for k in range(2): tmp[k]["p_holm"] = ph[k]; tmp[k]["p_used"] = ph[k]
    w_ = LP[(LP.stim == stim) & (LP.version == "after Laplacian (CSD)") & (LP.region == "ROI")].pivot(index="sub", columns="eye", values="log10_power").dropna(); d_ = w_.EC - w_.EO
    tt_ = _st.ttest_rel(w_.EC, w_.EO)
    tmp.append(dict(stimulation=stim, reference="Surface Laplacian", roi="sensorimotor", n=len(d_), t=tt_.statistic, df=len(d_) - 1, dz=d_.mean() / d_.std(ddof=1), p_uncorrected=tt_.pvalue, p_holm=np.nan, p_used=tt_.pvalue))
    for r_ in tmp:
        st_ = stars(r_["p_used"]); r_["label_text"] = ("Occipital" if r_["roi"] == "occipital" else "Sensorimotor") + f" dz = {r_['dz']:.2f}" + (" n.s." if st_ == "n.s." else st_)
        STAT[(stim, "average reference" if r_["reference"] == "Average reference" else "Laplacian", r_["roi"])] = r_; ROWS6.append(r_)
pd.DataFrame(ROWS6)[["stimulation", "reference", "roi", "n", "t", "df", "dz", "p_uncorrected", "p_holm", "label_text"]].to_csv(_A3 / "Figure6A_ROI_stats.csv", index=False, encoding="utf-8-sig")
def roi_geom(grp, m=0.05):     # same ellipse as roi_outline()
    g = np.array([POS[l] for l in grp]); c = (g.max(0) + g.min(0)) / 2; w, h = g.max(0) - g.min(0) + 2 * m; return c, w * 1.25, h * 1.25 + 0.02
AXA = {}
for i, stim in enumerate(("Real", "Sham")):
    for j, ver in enumerate(("average reference", "Laplacian")):
        ax = fig.add_axes([0.03 + j * 0.25, 0.70 - i * 0.215, 0.23, 0.20]); g = T[(T.stim == stim) & (T.version == ver)]
        im = topo(ax, list(g.channel), g.t.values, vmax, sig=None, rois=False); AXA[(stim, ver)] = ax              # rev1108: no FDR circles
        for reg_, grp_ in ((("occipital", OCC), ("sensorimotor", SMR)) if ver == "average reference" else (("sensorimotor", SMR),)):
            c_, W_, H_ = roi_geom(grp_); sig_ = STAT[(stim, ver, reg_)]["p_used"] < .05
            ax.add_patch(Ellipse(c_, W_, H_, fill=False, ls="-" if sig_ else "--", lw=1.5 if sig_ else .8, color="k", gid=f"roi_{reg_}"))
        ax.set_title(("Average reference" if j == 0 else "Surface Laplacian") if i == 0 else "", fontsize=9)
        if j == 0: ax.text(-0.08, 0.5, f"{stim} TMS", transform=ax.transAxes, rotation=90, va="center", ha="center", fontsize=9)
        if i == 0 and j == 0: plabel(ax, "A", dx=-0.05)
cb = fig.add_axes([0.53, 0.53, 0.010, 0.30]); fig.colorbar(im, cax=cb); cb.set_title("t (EC − EO)", fontsize=8); cb.tick_params(labelsize=8)
panels = (("occ", "Occipital" + chr(10) + "(average reference)"), ("smr", "Sensorimotor" + chr(10) + "(Laplacian)"))
PB = []
for i, stim in enumerate(("Real", "Sham")):
    for j, (reg, lab) in enumerate(panels):
        ax = fig.add_axes([0.66 + j * 0.18, 0.745 - i * 0.215, 0.13, 0.16])
        if reg == "occ":   # round 11: average reference, Holm over the two regions (manipulation check)
            w_ = AV[(AV.stim == stim) & (AV.region == "occipital")].pivot(index="sub", columns="eye", values="log10_power").dropna()
            r_ = MAN[(MAN.check == "3a") & (MAN.stim == stim) & (MAN.method == "Welch") & (MAN.region == "occipital")].iloc[0]; p = r_.p_holm_2
        else:              # unchanged: after the Laplacian, uncorrected (volume-conduction check)
            w_ = LP[(LP.stim == stim) & (LP.version == "after Laplacian (CSD)") & (LP.region == "ROI")].pivot(index="sub", columns="eye", values="log10_power").dropna()
            r_ = MAN[(MAN.check == "3b") & (MAN.stim == stim) & (MAN.method == "Welch, after Laplacian (CSD)") & (MAN.region == "ROI")].iloc[0]; p = r_.p
        PB.append(dict(stim=stim, panel=reg, n=len(w_), dz=r_.dz, p_plotted=p, stars=stars(p)))
        for _, r in w_.iterrows(): ax.plot([0, 1], [r.EO, r.EC], color="0.75", lw=.6, marker="o", ms=2, mfc="0.5", mec="none")
        ax.plot([0, 1], [w_.EO.mean(), w_.EC.mean()], color="k", lw=1.6, marker="o", ms=4)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["EO", "EC"] if i == 1 else ["", ""], fontsize=8); ax.set_xlim(-.35, 1.35); ax.tick_params(labelsize=7)
        ax.text(.5, .98, stars(p), transform=ax.transAxes, ha="center", va="top", fontsize=9)
        if i == 0: ax.set_title(lab, fontsize=8)
        if j == 0: ax.set_ylabel(stim + " TMS" + chr(10) + "log10 power (µV²/Hz)", fontsize=7.5)
        else: ax.set_ylabel("log10 Laplacian power (a.u.)", fontsize=7.5, labelpad=1)
        if i == 0 and j == 0: plabel(ax, "B", dx=-0.45)
def dzplot(ax, measures, title, letter):
    xs, labs = [], []
    for k, ms in enumerate(measures):
        for li, lv in enumerate(io.LEV):
            r = TO[(TO.measure == ms) & (TO.intensity == lv)].iloc[0]; x = k * 3.6 + li
            bd = TB5[(TB5.measure == ms) & (TB5.intensity == lv)].smallest_passing_bound_dz.iloc[0]
            ax.add_patch(plt.Rectangle((x - .38, -bd), .76, 2 * bd, color="0.88", lw=0, zorder=0))
            dz_, lo_, hi_ = -r.dz, -r.ci90_high / r.sd_diff, -r.ci90_low / r.sd_diff          # round 15: EC − EO
            ax.errorbar(x, dz_, yerr=[[dz_ - lo_], [hi_ - dz_]], fmt="o", color=CL3[li], ms=4, capsize=2, zorder=3)
        xs.append(k * 3.6 + 1); labs.append(ms.replace("ERM ", "").replace("ITPC ", ""))
    ax.axhline(0, color=GREY, lw=.6); ax.set_xticks(xs); ax.set_xticklabels(labs, fontsize=8); ax.set_ylim(-1.3, 1.3); ax.set_xlim(-.8, (len(measures) - 1) * 3.6 + 2.8)
    ax.set_ylabel(f"EC − EO ({_DZ}, 90% CI)", fontsize=8); ax.set_title(title, fontsize=9); plabel(ax, letter)
axC = fig.add_axes([0.07, 0.07, 0.15, 0.34]); dzplot(axC, ["log MEP"], "MEP", "C")
axD = fig.add_axes([0.31, 0.07, 0.38, 0.34]); dzplot(axD, ["ERM P30", "ERM N45", "ERM P60", "ERM N100"], "TEP components (ERM)", "D")
axE = fig.add_axes([0.78, 0.07, 0.20, 0.34]); dzplot(axE, ["ITPC IMF6", "ITPC IMF7"], "ITPC", "E")
for lv, col in zip(io.LEV, CL3): axE.plot([], [], "o", color=col, label=f"{lv}% RMT")
axE.add_patch(plt.Rectangle((0, 0), 0, 0, color="0.88", label="smallest TOST bound")); axE.legend(frameon=False, fontsize=6.5, loc="lower left")
# ---- rev1108: ROI labels outside the head, two lines, 6.5 pt, leader line (0.4 pt) to the outline (sensorimotor: upper left; occipital: lower right)
fig.canvas.draw(); _rr = fig.canvas.get_renderer(); LABS6 = []
for (stim, ver), ax in AXA.items():
    Lm = ax.get_xlim()[1]; rH = Lm / 1.02; bb = ax.get_window_extent(_rr); din = 2 * Lm / (bb.height / fig.dpi); mar = .03 * din; lh = 6.5 * 1.2 / 72 * din
    hw = lambda y: np.sqrt(max(rH ** 2 - y ** 2, 0))
    for reg_, grp_ in ((("sensorimotor", SMR), ("occipital", OCC)) if ver == "average reference" else (("sensorimotor", SMR),)):
        r_ = STAT[(stim, ver, reg_)]; l1, l2 = r_["label_text"].split(" dz = "); l2 = _DZ + " = " + l2   # 1114
        t1 = ax.text(0, 0, l1, fontsize=6.5, clip_on=False, gid=f"roilab_{reg_}"); t2 = ax.text(0, 0, l2, fontsize=6.5, clip_on=False, gid=f"roilab_{reg_}")
        w1 = t1.get_window_extent(_rr).width / fig.dpi * din; w2 = t2.get_window_extent(_rr).width / fig.dpi * din
        c_, W_, H_ = roi_geom(grp_)
        if reg_ == "sensorimotor":   # upper left, left-aligned block
            ytop = Lm - .02 * din; y1b = ytop - lh; y2b = y1b - lh
            xl = min(-hw(y1b) - mar - w1, -hw(y2b) - mar - w2)
            t1.set_position((xl, ytop)); t1.set_va("top"); t1.set_ha("left"); t2.set_position((xl, y1b)); t2.set_va("top"); t2.set_ha("left")
            p0 = (xl + w2 + .01 * din, y1b - lh / 2); a_ = np.deg2rad(135)
        else:                        # lower right, right-aligned block
            ybot = -Lm + .02 * din; y2t = ybot + lh; y1t = y2t + lh
            xr = max(hw(y1t) + mar + w1, hw(y2t) + mar + w2)
            t1.set_position((xr, y2t)); t1.set_va("bottom"); t1.set_ha("right"); t2.set_position((xr, ybot)); t2.set_va("bottom"); t2.set_ha("right")
            p0 = (xr - w1 - .01 * din, y2t + lh / 2); a_ = np.deg2rad(-45)
        p1 = (c_[0] + W_ / 2 * np.cos(a_), c_[1] + H_ / 2 * np.sin(a_))
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="k", lw=.4, clip_on=False, gid=f"leader_{reg_}")
        LABS6.append(dict(stimulation=stim, reference=ver, roi=reg_, text=r_["label_text"], ellipse_solid=bool(r_["p_used"] < .05)))
# rev1108 output: 600-dpi render with the original call (checks), 180 mm PNG, PDF; geometry of panels B-E for the pixel comparison
import json as _js
from PIL import Image as _Im
fig.savefig(_A3 / "_check" / "Figure6_rev_600dpi_tight.png", dpi=600, bbox_inches="tight")
_r = fig.canvas.get_renderer(); _bb = fig.get_tightbbox(_r); _x0, _y1 = _bb.x0 - .1, _bb.y1 + .1
def _px(b):   # display bbox (pixels at fig.dpi) -> pixel box of the 600-dpi tight render (x0, y0 from the top, x1, y1)
    return [int(np.floor((b.x0 / fig.dpi - _x0) * 600)), int(np.floor((_y1 - b.y1 / fig.dpi) * 600)), int(np.ceil((b.x1 / fig.dpi - _x0) * 600)), int(np.ceil((_y1 - b.y0 / fig.dpi) * 600))]
GEO = dict(A=[_px(a.get_tightbbox(_r)) for a in AXA.values()], BE=[_px(a.get_tightbbox(_r)) for a in fig.axes if a not in AXA.values() and (a.get_position().x0 > .6 or a.get_position().y1 < .45)],
           colorbar=[_px(a.get_tightbbox(_r)) for a in fig.axes if a not in AXA.values() and .5 < a.get_position().x0 < .6], labels=LABS6,
           label_boxes=[dict(stim=s_, ver=v_, gid=t.get_gid(), text=t.get_text(), box=_px(t.get_window_extent(_r))) for (s_, v_), a in AXA.items() for t in a.texts if (t.get_gid() or "").startswith("roilab")])
_js.dump(GEO, open(_A3 / "_check" / "fig6_geometry.json", "w"), indent=1)
_win = _bb.width + 2 * .1
_BBJ = _P(data_dir + r"\Revision_2026\Manuscript_Figures\Figures_2026-11-14_minor") / "_check" / "fig6_tightbbox.json"   # 1114
import io as _bio, matplotlib.figure as _mf; from matplotlib.transforms import Bbox as _Bb; _gtb = _mf.Figure.get_tightbbox; _REC = []   # 1114
def _rec(self, *a, **k): b = _gtb(self, *a, **k); _REC.append(b); return b   # 1114
_mf.Figure.get_tightbbox = _rec
if M1114 == "rev": _o = _js.load(open(_BBJ)); _win = _o["win"]; _FIXBB = _Bb.from_extents(*_o["save_tight"]).padded(.1)   # 1114: original canvas
fig.savefig(_bio.BytesIO(), format="png", dpi=4252 / _win, bbox_inches="tight", facecolor="white"); _mf.Figure.get_tightbbox = _gtb   # 1114: probe
_tb = [float(v) for v in _REC[0].extents]
if M1114 == "repro": _js.dump(dict(win=_win, save_tight=_tb), open(_BBJ, "w"), indent=1)   # 1114
else:
    _js.dump(dict(original=_o, new_save_tight=_tb), open(_BBJ.with_name("fig6_tightbbox_rev.json"), "w"), indent=1)
    assert _tb[0] > _FIXBB.x0 and _tb[1] > _FIXBB.y0 and _tb[2] < _FIXBB.x1 and _tb[3] < _FIXBB.y1, "new content outside the original canvas"
_ST = _P(data_dir + r"\Revision_2026\Manuscript_Figures\Figures_2026-11-14_minor") / ("_check/repro/repro_Figure6_rev" if M1114 == "repro" else "Figure6_1114")   # 1114
fig.savefig(f"{_ST}.png", dpi=4252 / _win, bbox_inches="tight" if M1114 == "repro" else _FIXBB, facecolor="white"); fig.savefig(f"{_ST}.pdf", bbox_inches="tight" if M1114 == "repro" else _FIXBB, facecolor="white")
_im = _Im.open(f"{_ST}.png").convert("RGB")
if _im.size[0] != 4252: _im = _im.resize((4252, round(_im.size[1] * 4252 / _im.size[0])), _Im.LANCZOS)
_im.save(f"{_ST}.png", dpi=(600, 600)); plt.close(fig)   # 1114
_FIG6_AXA = AXA

# ---- rev1108: Supplementary figure = the original panel A (per-channel FDR circles, dashed ROI outlines), same colour scale
figS = plt.figure(figsize=(4.6, 4.3))
figS.suptitle("Per-channel EC − EO t values of pre-stimulus alpha power" + chr(10) + "(FDR q < .05 over 61 channels)", fontsize=10, fontweight="bold", y=.995)
for i, stim in enumerate(("Real", "Sham")):
    for j, ver in enumerate(("average reference", "Laplacian")):
        ax = figS.add_axes([0.07 + j * 0.40, 0.47 - i * 0.44, 0.38, 0.38]); g = T[(T.stim == stim) & (T.version == ver)]
        im = topo(ax, list(g.channel), g.t.values, vmax, sig=(g.q_fdr61 < .05).values)
        if i == 0: ax.set_title("Average reference" if j == 0 else "Surface Laplacian", fontsize=9)
        if j == 0: ax.text(-0.08, 0.5, f"{stim} TMS", transform=ax.transAxes, rotation=90, va="center", ha="center", fontsize=9)
cbS = figS.add_axes([0.88, 0.22, 0.018, 0.50]); figS.colorbar(im, cax=cbS); cbS.set_title("t (EC − EO)", fontsize=8); cbS.tick_params(labelsize=8)
figS.savefig(_A3 / "FigureS_alpha_perchannel.png", dpi=600, bbox_inches="tight", facecolor="white"); figS.savefig(_A3 / "FigureS_alpha_perchannel.pdf", bbox_inches="tight", facecolor="white")
_im = _Im.open(_A3 / "FigureS_alpha_perchannel.png").convert("RGB"); _im.save(_A3 / "FigureS_alpha_perchannel.png", dpi=(600, 600)); plt.close(figS)

print(pd.DataFrame(PB).to_string(index=False))
