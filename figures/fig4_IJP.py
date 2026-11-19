data_dir = r"path/to/data"   # release: set this path
# _IJP copy (user 2026-11-10, Part F2) of fig4_1111.py made by Figures_IJP_2026-11-10/code/make_IJP_copies.py: figure title 'Early cortical responses ...' -> 'Early EEG responses ...'
# (IJP_MODE=ijp, default); IJP_MODE=repro = Figure4_1111 unchanged. Output: Figures_IJP_2026-11-10/Figure4_IJP_rgb.png/.pdf (repro: _check/repro_Figure4_1111). Lines marked # IJP.
# _1111 copy (user 2026-10-28, Figure4_2026-11-11) made by make_fig4_1111.py from Figures_3_4_6_2026-11-08/code/fig4_p6_rev1108.py: panel B legend with Line2D handles (F4_MODE=rev); F4_MODE=repro = unchanged drawing.
# _rev1108 copy (user 2026-10-27, task 2) of qc/figures_p6/fig4_p6.py (the original; the rev1107 error bars, y limit and legend are not used),
# made by Figures_3_4_6_2026-11-08/code/make_rev1108_copies.py. Drawing unchanged; only the output: Figure4_rev.png (180 mm, 600 dpi, RGB), .pdf, and a
# 600-dpi render with the original settings in _check/ for the pixel comparison with Fig4_p6.png.
# Part 6 copy made by qc/figures_p6/make_p6_copies.py (2026-10-23): k/n labels above the bars removed; nothing else changed.
# A3 copy made by qc/figures_a3/make_a3_figure_copies.py (part 4, 2026-10-22).
import sys as _sys; _sys.path.insert(0, data_dir + r'\Revision_2026\Manuscript_Figures\code')
from pathlib import Path as _P
"""Round 15 (user prompt 2026-10-27, section 2): Figure 4 (R2-1) redrawn from data/R21_r15_cells.csv (r21_round15.py).
Layout of Fig8_round8 (fig_round8.py) kept; changes: (A) k/n above every bar, zero-valued bars get a short mark on the baseline so all four
conditions are visible at every intensity; (B) log SNR with the local-baseline null.
VERSION=local (default; main) or orig (previous null, for the comparison in the report/explainer only).
Outputs: Fig4_r15.png/.tiff (local), Fig4_r15_orig_null.png (orig)."""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import tep_io as io

V = os.environ.get("VERSION", "local"); D, MF = io.DATA, io.ROOT
_IJP = os.environ.get("IJP_MODE", "ijp"); assert _IJP in ("repro", "ijp")   # IJP
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 9, "axes.linewidth": .8})
CEO, CEC, GREY = "#E08A2C", "#3B8C8C", "0.55"
def plabel(ax, s): ax.text(-0.12, 1.06, s, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
C = pd.read_csv(D / "R21_r15_cells.csv"); C = C[C["sub"].isin(io.MAIN17)]
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.3)); fig.suptitle(("Early cortical responses with real and sham TMS at three intensities" if _IJP == "repro" else "Early EEG responses with real and sham TMS at three intensities"),   # IJP
    fontsize=11, fontweight="bold", y=1.03)   # round 15: previous title no longer matches the local-null counts [TO CONFIRM by the user]
w = .18
for k, (stim, eye, col, hatch) in enumerate((("Real", "EO", CEO, None), ("Real", "EC", CEC, None), ("Sham", "EO", CEO, "////"), ("Sham", "EC", CEC, "////"))):
    g = C[(C.stim == stim) & (C.eye == eye)].groupby("intensity")[f"responder_{V}"].agg(["sum", "count"]).sort_index()
    xs = np.arange(3) + (k - 1.5) * w; val = g["sum"] / g["count"]
    ax[0].bar(xs, val, w, color=col if stim == "Real" else "white", edgecolor=col, hatch=hatch, lw=1, label=f"{stim} {eye}")
    for x_, v_, kk, nn in zip(xs, val, g["sum"], g["count"]):
        if kk == 0: ax[0].plot([x_ - w / 2 + .01, x_ + w / 2 - .01], [.006, .006], color=col, lw=2.2, solid_capstyle="butt", clip_on=False)
ax[0].set_xticks(range(3)); ax[0].set_xticklabels(["80%", "100%", "120%"]); ax[0].set_xlabel("Intensity (% RMT)"); ax[0].set_ylabel("Proportion with a reliable response")
ax[0].set_ylim(0, 1.45); ax[0].set_yticks(np.arange(0, 1.01, .2))
ax[0].legend(frameon=False, fontsize=7.5, ncol=4, loc="upper center", handlelength=1.2, columnspacing=.8); plabel(ax[0], "A")
for stim, ls in (("Real", "-"), ("Sham", "--")):
    for eye, col in (("EO", CEO), ("EC", CEC)):
        g = C[(C.stim == stim) & (C.eye == eye)].groupby("intensity")[f"log_snr_{V}"].agg(["mean", "std", "count"])
        ax[1].errorbar(g.index + (-3 if eye == "EO" else 3), g["mean"], yerr=g["std"] / np.sqrt(g["count"]), fmt="o", ls=ls, color=col, ms=4, capsize=2, label=f"{stim} {eye}")
ax[1].axhline(0, color=GREY, lw=.6); ax[1].set_xticks([80, 100, 120]); ax[1].set_xlabel("Intensity (% RMT)"); ax[1].set_ylabel("log SNR (15–60 ms GMFP)")
F4 = os.environ.get("F4_MODE", "rev"); assert F4 == "rev"   # IJP: legend of Figure4_1111
if F4 == "repro": ax[1].legend(frameon=False, fontsize=7.5)
else:   # 1111: Line2D handles with the line style of the curves (errorbar handles did not show solid vs dashed)
    from matplotlib.lines import Line2D
    _ln = {c.get_label(): c.lines[0] for c in ax[1].containers}
    _h = [Line2D([], [], color=_ln[l].get_color(), ls=_ln[l].get_linestyle(), lw=_ln[l].get_linewidth(), marker=_ln[l].get_marker(), ms=_ln[l].get_markersize()) for l in _ln]
    # handlelength 3.0 (default 2.0): at 2.0 the marker covers the only gap of the dash pattern and the sham handle looks solid; loc fixed to
    # "center right" = the location that loc="best" chose in the original (legend right edge and vertical centre unchanged)
    ax[1].legend(_h, list(_ln), frameon=False, fontsize=7.5, loc="center right", handlelength=3.0)
plabel(ax[1], "B")
fig.tight_layout()
stem = "Figure4_rev"; assert V == "local"
from PIL import Image as _Im
_OUTD = _P(data_dir + r"\Revision_2026\Manuscript_Figures\Figures_IJP_2026-11-10"); (_OUTD / "_check").mkdir(parents=True, exist_ok=True)   # IJP
if _IJP == "repro": _OUTD = _OUTD / "_check"   # IJP
stem = "Figure4_IJP_rgb" if _IJP == "ijp" else "repro_Figure4_1111"   # IJP
print("legend B bbox (px at fig dpi):", [round(v, 2) for v in ax[1].get_legend().get_window_extent(fig.canvas.get_renderer()).extents])
_r = fig.canvas.get_renderer(); _bb = fig.get_tightbbox(_r); _win = _bb.width + 2 * .1                                  # bbox_inches='tight' adds pad_inches = 0.1
fig.savefig(_OUTD / f"{stem}.png", dpi=4252 / _win, bbox_inches="tight", facecolor="white"); fig.savefig(_OUTD / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
_im = _Im.open(_OUTD / f"{stem}.png").convert("RGB")
if _im.size[0] != 4252: _im = _im.resize((4252, round(_im.size[1] * 4252 / _im.size[0])), _Im.LANCZOS)
_im.save(_OUTD / f"{stem}.png", dpi=(600, 600)); plt.close(fig); print("written", stem, _im.size)
