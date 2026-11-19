data_dir = r"path/to/data"   # release: set this path
# _1114 copy (user 2026-11-14, Figures_2026-11-14_minor) of Review_2026-11-10/code/fig3_1110.py, always drawn as B11_MODE=rev (= Figure3_1110).
# M1114=repro: unchanged -> _check/repro; M1114=rev: panel B title 'ERM (CEEMDAN IMFs 4–7)'. Panel A untouched. Lines marked # 1114.
# _1110 copy (user 2026-10-28, Review_2026-11-10, B11) made by make_fig3_1110.py from Figures_3_4_6_2026-11-08/code/fig3_s9_s10_c4_rev1108.py: 0-10 ms black bar in panel B (B11_MODE=rev); B11_MODE=repro = unchanged drawing.
# _rev1108 copy (user 2026-10-27, task 2) of Figures_3_4A_2026-11-07/code/fig3_s9_s10_c4_rev1107.py, made by Figures_3_4_6_2026-11-08/code/make_rev1108_copies.py.
# Changes (drawing only): brackets removed; component labels (text only) at the peaks: placed at the time of the maximum (P30 / P60) or minimum
# (N45 / N100) of the six curves within the component window of the panel (A: conventional TEP windows, B: ERM windows, data/TEP_component_windows.csv),
# 4 pt above / below the extreme of the curves under the label; moved within the window when needed to avoid the axes frame, the red bars or the other
# label (P60 vs P30, 3 pt). Everything else as rev1107. Output Figures_3_4_6_2026-11-08.
# _rev1107 copy (user 2026-10-27), made by Figures_3_4A_2026-11-07/code/make_rev1107_copies.py: Figure 3 only. Data, statistics, colours, axis ranges
# and axes layout unchanged. Changes (drawing only): 1. grey component windows (axvspan 0.88) removed; 2. component windows marked by thin black
# brackets (1 pt, short end ticks) at the top of each panel, labels above them (panel B keeps the stars); N45 on an upper tier, P30 / P60 / N100 on a
# lower tier with P60 offset 4 pt downwards, P30 / P60 labels moved 3 pt outwards (overlaps checked in check_figures_rev1107.py); 3. black removal window (-5..10 ms) kept; 4. legend labels
# with spaces; 5. output Figures_3_4A_2026-11-07/Figure3_rev.png (180 mm, 600 dpi, RGB) and .pdf. FigS10 / FigS9 of the original are not produced.
# C4 copy made by qc/figures_c4/make_c4_copies.py (part 4, 2026-10-23); outputs in qc/figures_c4/.
import sys as _sys; _sys.path.insert(0, data_dir + r'\Revision_2026\Manuscript_Figures\code')
from pathlib import Path as _P
import os as _os
_MODE = 'rev'   # 1114: always the Figure3_1110 drawing
M1114 = _os.environ.get('M1114', 'rev'); assert M1114 in ('repro', 'rev')   # 1114
_OUTD = _P(data_dir + r'\Revision_2026\Manuscript_Figures\Figures_2026-11-14_minor') / ('_check/repro' if M1114 == 'repro' else '')   # 1114
_CKD = _P(data_dir + r'\Revision_2026\Manuscript_Figures\Figures_2026-11-14_minor') / '_check' / 'fig3'   # 1114
"""Fig 5 (master plan C; decisions 2026-10-16), format of the original Fig 5 (spTMS_ERM.m plot block: six condition colours, 3-pt lines,
black removal window, Times New Roman bold, legend outside, ROI head inset), two rows, 17 participants:
  A  conventional TEP (_F_TEP, 0.1-30 Hz), -100..400 ms
  B  ERM (CEEMDAN IMF4-7), 0..400 ms, red bars = 20-ms windows with an FDR-significant intensity main effect (sliding-window rmANOVA,
     data/TEP_slidingwindow_anova_ERM.csv)
Both rows: P30/N45/P60 labels at the grand-average peaks and light-grey windows (peak +-10 ms; data/TEP_component_windows.csv); N100 window 120-140 ms.
Also FigS_Fig5B_uncorrected (uncorrected bars, Supplementary) and FigS_SE_slidingwindow (S-E: partial eta squared of the intensity effect per window
for ERM, band-pass control and conventional TEP, FDR windows marked).
Outputs (600 dpi PNG + TIFF): Fig5_new, FigS_Fig5B_uncorrected, FigS_SE_slidingwindow; data/Fig5_plot_data.csv
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
import tep_io as io

S17 = io.MAIN17; COLS = np.array([[249, 149, 78], [249, 114, 19], [197, 85, 5], [122, 184, 184], [109, 173, 157], [57, 106, 106]]) / 255
CONDS = [(e, lv) for e in ("EO", "EC") for lv in io.LEV]
W = pd.read_csv(io.DATA / "TEP_roi_waveforms.csv"); W = W[W.file == "TEP"]; tc = [c for c in W.columns if c.startswith("t_")]; t_tep = np.arange(-2000, -2000 + len(tc))
TEP = {(int(r.id[6:9]), r.id[-2:], int(r.intensity)): np.asarray(r[tc], float) for _, r in W.iterrows() if r.id.split("_")[2] == "Real"}
ERM = {(s, e): io.load_erm(s, e) for s in S17 for e in io.EYES}; t_erm = np.round(ERM[(S17[0], "EO")]["time"]).astype(int)
GA = {"TEP": {c: np.mean([TEP[(s, c[0], c[1])] for s in S17], 0) for c in CONDS},
      "ERM": {c: np.mean([ERM[(s, c[0])]["erm_roi"][io.LEV.index(c[1])] for s in S17], 0) for c in CONDS}}
CW = pd.read_csv(io.DATA / "TEP_component_windows.csv")
# 2026-10-21 (round 8): stars of the intensity main effect of the ERM component rmANOVA next to the labels of panel B (GG p; * < .05, ** < .01, *** < .001)
_CA = pd.read_csv(io.DATA / "TEP_component_anova.csv"); _CA = _CA[_CA.data == "ERM"].set_index("component").p_intensity
STAR = {c: ("***" if p < .001 else "**" if p < .01 else "*" if p < .05 else "") for c, p in _CA.items()}
TITLE = "TMS-evoked responses increase with stimulation intensity in both eye states"
SW = {k: pd.read_csv(io.DATA / f"TEP_slidingwindow_anova_{k}.csv") for k in ("ERM", "TEP", "BP")}
rows = [dict(data=k, eye=c[0], intensity=c[1], time_ms=int(tt), uV=float(v)) for k, tvec in (("TEP", t_tep), ("ERM", t_erm)) for c in CONDS
        for tt, v in zip(tvec, GA[k][c]) if -100 <= tt <= 400]
pass   # C4 copy: data/Fig5_plot_data.csv not rewritten

plt.rcParams.update({"font.family": "Times New Roman", "font.weight": "bold", "axes.labelweight": "bold", "axes.linewidth": 1.2})
ROIPOS = {"FC3": (-60.5, .227), "FC1": (-43.3, .141), "C3": (-89.3, .202), "C1": (-89.3, .100)}      # theta (deg), radius (waveguard64_woRef.ced)
def head(ax):
    ax.add_patch(Circle((0, 0), .5, fill=False, lw=1)); ax.add_patch(Polygon([[-.05, .495], [0, .56], [.05, .495]], closed=False, fill=False, lw=1))
    for th, r in ROIPOS.values(): ax.plot(r * np.sin(np.deg2rad(th)), r * np.cos(np.deg2rad(th)), "ko", ms=4)
    ax.set_xlim(-.6, .6); ax.set_ylim(-.6, .62); ax.set_aspect("equal"); ax.axis("off")
def panel(ax, kind, xlim, bars, label):
    tvec = t_tep if kind == "TEP" else t_erm; w = (tvec >= xlim[0]) & (tvec <= xlim[1])
    allv = np.concatenate([GA[kind][c][w] for c in CONDS]); lo, hi = np.floor(allv.min() - 1), np.ceil(allv.max() + 2.5)
    hs = [ax.plot(tvec[w], GA[kind][c][w], color=COLS[i], lw=1.5, zorder=2)[0] for i, c in enumerate(CONDS)]
    if xlim[0] < 0: ax.axvspan(-5, 10, color="k", lw=0, zorder=3)
    elif _MODE == "rev": ax.axvspan(max(-5, xlim[0]), 10, color="k", lw=0, zorder=3)   # 1110 (B11): 0-10 ms bar in panel B, style of panel A
    # rev1108: labels at the peaks (no brackets): extreme of the six curves within the component window, label 4 pt above (P) / below (N)
    PTU = (hi - lo) / (ax.get_position().height * ax.figure.get_figheight() * 72)          # data units per point
    for r in CW[CW.data == kind].itertuples():
        a0, a1 = r.win_start_ms, r.win_end_ms; ww = (tvec >= a0) & (tvec <= a1); M_ = np.array([GA[kind][c][ww] for c in CONDS]); tt_ = tvec[ww]
        pos_ = r.component.startswith("P"); k_ = np.unravel_index(np.argmax(M_) if pos_ else np.argmin(M_), M_.shape); x0_, y0_ = float(tt_[k_[1]]), float(M_[k_])
        ax.text(x0_, y0_ + (4 if pos_ else -4) * PTU, r.component + (STAR.get(r.component, "") if kind == "ERM" else ""), ha="center", va="bottom" if pos_ else "top",
                fontsize=12, fontweight="bold", zorder=4, gid=f"lab_{r.component}")
        LABPOS.append(dict(panel="A" if kind == "TEP" else "B", data=kind, component=r.component, window_ms=f"{a0}-{a1}", peak_time_ms=x0_, peak_uV=y0_,
                           peak_condition=f"{CONDS[k_[0]][1]}% RMT ({CONDS[k_[0]][0]})", label_anchor_uV=y0_ + (4 if pos_ else -4) * PTU, dx_pt=0.0))
    hb = None
    for a0 in bars: hb = ax.plot([a0, a0 + 20], [hi - .4, hi - .4], color=np.array([173, 52, 62]) / 255, lw=4, solid_capstyle="butt", zorder=4)[0]
    ax.set_xlim(*xlim); ax.set_ylim(lo, hi); ax.set_ylabel("Potential (µV)", fontsize=14); ax.tick_params(labelsize=11)
    ax.text(-0.07, 1.02, label, transform=ax.transAxes, fontsize=18, fontweight="bold")
    return hs, hb
def fig5(bar_col, stem, note):
    fig = plt.figure(figsize=(11, 9.0)); ax1 = fig.add_axes([0.08, 0.545, 0.66, 0.34]); ax2 = fig.add_axes([0.08, 0.065, 0.66, 0.34])
    fig.text(0.5, 0.975, TITLE, ha="center", va="top", fontsize=15, fontweight="bold")
    hs, _ = panel(ax1, "TEP", (-100, 400), [], "A")
    bars = SW["ERM"].win_start_ms[SW["ERM"][bar_col] < .05].values
    _, hb = panel(ax2, "ERM", (0, 400), bars, "B"); ax2.set_xlabel("Time (ms)", fontsize=14); ax1.set_xlabel("Time (ms)", fontsize=14)
    ax1.set_title("Conventional TEP (0.1–30 Hz)", fontsize=13); ax2.set_title("ERM (CEEMDAN IMF4–7)" if M1114 == "repro" else "ERM (CEEMDAN IMFs 4–7)", fontsize=13)   # 1114
    labs = ["80% RMT (EO)", "100% RMT (EO)", "120% RMT (EO)", "80% RMT (EC)", "100% RMT (EC)", "120% RMT (EC)"]
    hh = hs + ([hb] if hb is not None else []); ll = labs + (["TMS intensity effect"] if hb is not None else [])
    fig.legend(hh, ll, loc="center left", bbox_to_anchor=(0.76, 0.42), frameon=True, edgecolor="k", fontsize=11)
    head(fig.add_axes([0.80, 0.66, 0.14, 0.20]))
    # rev1108: final label placement. Each label starts at the time of the extreme within its window; its height is 4 pt beyond the extreme of the six
    # curves over the label's own horizontal extent (so the label never touches a curve). If the label would then touch the axes frame or a red bar, or
    # (P60) come closer than 3 pt to the P30 label, the label is moved in 1-ms steps within the component window, nearest positions first.
    import json as _js
    fig.canvas.draw(); _rr = fig.canvas.get_renderer()
    for _ax, _pn, _kind in ((ax1, "A", "TEP"), (ax2, "B", "ERM")):
        _tv = t_tep if _kind == "TEP" else t_erm; _G = np.array([GA[_kind][c] for c in CONDS]); _lo, _hi = _ax.get_ylim(); _x0l, _x1l = _ax.get_xlim()
        _ptu = (_hi - _lo) / (_ax.get_window_extent(_rr).height * 72 / fig.dpi)
        _red = [l.get_window_extent(_rr) for l in _ax.lines if l.get_linewidth() == 4]
        _L = {t.get_gid()[4:]: t for t in _ax.texts if (t.get_gid() or "").startswith("lab_")}; _placed = {}
        for _c in ("P30", "N45", "P60", "N100"):
            t = _L[_c]; d = [d for d in LABPOS if d["panel"] == _pn and d["component"] == _c][0]; a0, a1 = map(int, d["window_ms"].split("-")); pos_ = _c.startswith("P")
            wd = t.get_window_extent(_rr).width / _ax.get_window_extent(_rr).width * (_x1l - _x0l); ok = False
            _xy = [_ax.transData.transform(np.c_[np.interp(_xx := np.linspace(_tv[0], _tv[-1], 20 * len(_tv)), _tv, _tv), np.interp(_xx, _tv, _g)]) for _g in _G]
            def _dist(bb):   # shortest distance (pt) between the label box and the six curves, minus half the curve width (0.75 pt)
                return min(np.hypot(np.maximum(np.maximum(bb.x0 - q[:, 0], q[:, 0] - bb.x1), 0), np.maximum(np.maximum(bb.y0 - q[:, 1], q[:, 1] - bb.y1), 0)).min() for q in _xy) * 72 / fig.dpi - .75
            for extra in range(0, 9):          # extra distance (pt) beyond the 4 pt, only if no position in the window clears the curves by 3 pt
                for x_ in sorted(range(a0, a1 + 1), key=lambda v: (abs(v - d["peak_time_ms"]), v)):
                    span = (_tv >= x_ - wd / 2) & (_tv <= x_ + wd / 2); ext = _G[:, span].max() if pos_ else _G[:, span].min(); off = (4 + extra) * (1 if pos_ else -1)
                    t.set_position((x_, ext + off * _ptu)); bb = t.get_window_extent(_rr); ab = _ax.get_window_extent(_rr)
                    bad = bb.y0 < ab.y0 + 1 or bb.y1 > ab.y1 - 1 or bb.x0 < ab.x0 + 1 or bb.x1 > ab.x1 - 1 or _dist(bb) < 3
                    bad = bad or any(rb.x1 > bb.x0 and rb.x0 < bb.x1 and rb.y1 > bb.y0 - 2 and rb.y0 < bb.y1 + 2 for rb in _red)
                    bad = bad or any(min(bb.y1, o.y1) > max(bb.y0, o.y0) - 3 * fig.dpi / 72 and min(bb.x1, o.x1) > max(bb.x0, o.x0) - 3 * fig.dpi / 72 for o in _placed.values())
                    if not bad: ok = True; break
                if ok: break
            assert ok, (_pn, _c)
            _placed[_c] = t.get_window_extent(_rr); d.update(label_time_ms=float(x_), label_anchor_uV=float(ext + off * _ptu), extreme_over_label_uV=float(ext), offset_pt=abs(off))
    _js.dump(LABPOS, open(_CKD / f"fig3_label_positions_{_MODE}.json", "w"), indent=1)
    # rev1107: 180 mm wide at 600 dpi (4252 px; the figure is 11 in wide, so it is rendered at 4252 / 11 dpi and tagged 600 dpi), RGB, white; vector PDF
    from PIL import Image as _Im
    import pickle as _pk
    _D = _OUTD if _MODE == "rev" else _CKD   # 1110
    fig.savefig(_D / f"{stem}.png", dpi=4252 / 11, facecolor="white"); fig.savefig(_D / f"{stem}.pdf", facecolor="white")
    _im = _Im.open(_D / f"{stem}.png").convert("RGB"); assert _im.size[0] == 4252, _im.size; _im.save(_D / f"{stem}.png", dpi=(600, 600))
    _GEOM.append(fig)
_CKD.mkdir(parents=True, exist_ok=True); _GEOM = []; LABPOS = []
fig5("q_intensity", "repro_Figure3_1110" if M1114 == "repro" else "Figure3_1114", "FDR")   # 1114
print("done")
