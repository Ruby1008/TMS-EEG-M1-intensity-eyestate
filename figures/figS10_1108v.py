data_dir = r"path/to/data"   # release: set this path
# _1108v copy (user 2026-10-28, Figures_2026-11-08v) of Figures_2026-11-08t/code/figS10_1108t.py. Change (V_MODE=rev): the legend moves from the lower right of the 120% RMT panel (it covered participant 026) to the empty upper-left corner of the same panel; the script checks that the legend box overlaps no point, CI line or the zero line. V_MODE=repro reproduces FigureS10_1108t.png into _check/_scratch.
# _1108t copy (user 2026-10-28, Figures_2026-11-08t) of MEP_reprocess_2026-10-22/qc/single_trial_per_participant_plot_reprocessed.py (old S18, run with MEPSET=main_qc2_008). Changes (T_MODE=rev): BLUP diamonds and their legend entry removed (the plot script never fitted the model; it only read the BLUP column); output FigureS10_1108t.png/.pdf
# copy made by qc/make_b_copies.py (part 4, 2026-10-23); MEPSET = main_qc (Figure S18)
import sys as _sys; _sys.path.insert(0, data_dir + r'\Revision_2026\Manuscript_Figures\code')
"""Per-participant slopes of log MEP on the single-trial P30 (user 2026-10-20, task 1-(1)): OLS slope +- 95% CI per participant (MEP > 0 trials, EO and EC pooled),
BLUP slope of the random-slope model (open diamond), 100% and 120% RMT; filled circle = Spearman rho > 0. Times New Roman, 600 dpi, no n / correction wording.
Input: data/Fig6_single_trial_per_participant.csv. Output: FigS_single_trial_P30_MEP_per_participant.png/.tiff"""
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import tep_io as io
import os as _os
from pathlib import Path as _P
_SET = _os.environ.get("MEPSET", "main"); _RD = _P(data_dir + r"\Revision_2026\Manuscript_Figures\MEP_reprocess_2026-10-22") / _SET
_TM = _os.environ.get("T_MODE", "rev"); _OUT = _P(data_dir + r"\Revision_2026\Manuscript_Figures\Figures_2026-11-08v")   # 1108v
_VM = _os.environ.get("V_MODE", "rev")

P = pd.read_csv(_RD / "analysis" / "Fig6_single_trial_per_participant_reprocessed.csv")
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 9})
fig, ax = plt.subplots(1, 2, figsize=(7.0, 4.2), sharey=True)
order = P[P.intensity == 100].sort_values("ols_slope")["sub"].tolist()
for a, lv in zip(ax, (100, 120)):
    g = P[P.intensity == lv].set_index("sub").loc[order]; y = np.arange(len(g))
    ci = stats.t.ppf(.975, g.n_trials - 2) * g.ols_se
    a.axvline(0, color="0.5", lw=.7)
    a.errorbar(g.ols_slope, y, xerr=ci, fmt="none", ecolor="0.6", lw=.9)
    pos = g.rho > 0
    a.plot(g.ols_slope[pos], y[pos.values], "o", color="#2E6E6A", ms=4.5, label="Spearman ρ > 0")
    a.plot(g.ols_slope[~pos], y[~pos.values], "o", mfc="white", mec="#2E6E6A", ms=4.5, label="Spearman ρ ≤ 0")
    if _TM == "repro": a.plot(g.blup_slope, y, "D", mfc="none", mec="#C5550A", ms=3.5, label="BLUP slope")   # 1108t: random-slope BLUPs not shown
    a.set_yticks(y); a.set_yticklabels([f"{s:03d}" for s in g.index]); a.set_title(f"{lv}% RMT"); a.set_xlabel("Slope of log MEP on P30 (per µV)")
ax[0].set_ylabel("Participant"); leg = ax[1].legend(frameon=False, fontsize=7.5, loc="lower right" if _VM == "repro" else "upper left"); fig.tight_layout()
if _VM != "repro":   # 1108v: legend must not touch any point, CI line or the zero line of the 120% panel
    fig.canvas.draw(); r = fig.canvas.get_renderer(); lb = leg.get_window_extent(r); hits = []
    for art in ax[1].lines + ax[1].collections:      # point markers and the zero line (Line2D), CI lines (LineCollection)
        pad = (4.5 / 2 + 1) * fig.dpi / 72 if art in ax[1].lines and art.get_marker() == "o" else 1 * fig.dpi / 72
        segs = art.get_segments() if hasattr(art, "get_segments") else [np.column_stack(art.get_data())]
        for sg in segs:
            xy = art.get_transform().transform(np.asarray(sg, float))
            if art.get_linestyle() not in ("None", "") and len(xy) == 2: xy = np.linspace(xy[0], xy[1], 400)
            inside = (xy[:, 0] > lb.x0 - pad) & (xy[:, 0] < lb.x1 + pad) & (xy[:, 1] > lb.y0 - pad) & (xy[:, 1] < lb.y1 + pad)
            if inside.any(): hits.append(f"{type(art).__name__} {art.get_label()}")
    print("legend box (px):", [round(v) for v in lb.extents], "| overlaps:", hits or "none")
    print("axes boxes (px):", [[round(v) for v in a.get_window_extent(r).extents] for a in ax])
if _VM == "repro": fig.savefig(_OUT / "_check" / "_scratch" / "repro_FigureS10_1108t.png", dpi=600)
else: fig.savefig(_OUT / "FigureS10_1108v.png", dpi=600); fig.savefig(_OUT / "FigureS10_1108v.pdf")
print("saved")
