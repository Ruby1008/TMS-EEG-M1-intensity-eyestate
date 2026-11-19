# _1108t copy (user 2026-10-28, Figures_2026-11-08t) of Figures_2026-11-08n/code/figS3_oldS8_1108n.py (run with LBL_ALT=1, the final 'Potential' version). Changes (T_MODE=rev): figure title 'mean ± SEM' -> 'mean ± SE'; output FigureS3_1108t.png/.pdf
# _1108n copy (user 2026-10-28, Figures_2026-11-08n) of Figure_S2_2026-11-08/code/figS8_sd_single_subject_tep_rev1108c.py. Changes: alt only: y label 'µV' -> 'Potential (µV)' -> _check/alt_potential/FigureS3_1108n_alt.png (run with LBL_ALT=1).
# _rev1108c copy (user 2026-10-27, Figure_S2_2026-11-08) of Figures_all_2026-11-08/code/figS8_sd_single_subject_tep_rev1108.py: the semicolon of the title -> comma;
# outputs FigureS8_rev_v2.png/.pdf in Figure_S2_2026-11-08. Everything else unchanged.
# _rev1108 copy (user 2026-10-27, Figures_all_2026-11-08) of Manuscript_Figures/code/sd_single_subject_tep.py (FigS_SD_overview_EO120.png = Figure S8, MD5).
# Overview figure only (the per-participant PDF and the n-trials CSV are not written). Change (REV1108_MODE=rev): grey removal bar (0.85) of every panel ->
# black solid bar. Data, colours (blue SEM band), layout unchanged. Output FigureS8_rev.png/.pdf; repro: _check/repro_FigS8.png.
"""S-D single-participant TEPs (master plan E; user 2026-10-18, section 4.2).
Data: Preproc\\Final\\<id>_F_ERM.set (0.1-80 Hz, before the 30 Hz low-pass; average reference; baseline -100..-6 ms), real stimulation.
One page per participant: rows EO / EC; columns 80, 100, 120% RMT = ROI (FC3/FC1/C3/C1) mean ± SEM across trials, and a butterfly plot at 120% RMT
(61 scalp channels, ROI mean in black); -100..400 ms; grey: removed and interpolated window (-5 to 10 ms).
Pages: the 17 main participants first, then the other participants with at least one available real-stimulation condition (marked).
Overview: all available participants, real EO 120% RMT, ROI mean ± SEM.
Outputs: Manuscript_Figures\\FigS_SD_single_participant_TEPs.pdf, FigS_SD_overview_EO120.png/.tiff (600 dpi), data/SD_single_participant_ntrials.csv"""
data_dir = r"path/to/data"   # release: set this path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys as _s0; _s0.path.insert(0, data_dir + r'\Revision_2026\Manuscript_Figures\code')
import tep_io as io
import sys as _s; _s.path.insert(0, data_dir + r'\Revision_2026\Manuscript_Figures\Figures_all_2026-11-08\code'); import rev1108_helpers as H
import os as _os; from pathlib import Path
CMODE = _os.environ.get("REV1108C_MODE", "rev"); assert CMODE in ("rev", "repro")   # rev1108c
assert H.MODE == "rev"   # the rev1108 drawing is the starting point

plt.rcParams.update({"font.family": "Times New Roman", "font.size": 7}); COL = {80: "#9ecae1", 100: "#4292c6", 120: "#08306b"}
others = sorted((set(io.AVAIL["EO"]) | set(io.AVAIL["EC"])) - set(io.MAIN17)); order = [(s, True) for s in io.MAIN17] + [(s, False) for s in others]
def load(s, e):
    if s not in io.AVAIL[e]: return None
    x, labs, t, inten, _ = io.load_set(io.FINAL, f"spTMS_{s:03d}_Real_{e}_F_ERM"); t = np.round(t).astype(int)
    sc = [i for i, l in enumerate(labs) if l not in ("M1", "M2")]; ri = [labs.index(c) for c in io.ROI]; w = (t >= -100) & (t <= 400)
    return dict(t=t[w], roi=x[ri][:, w].astype(float).mean(0), all=x[sc][:, w].astype(float), inten=inten)
rows, ov = [], {}
# rev1108: the per-participant PDF pages and data/SD_single_participant_ntrials.csv are not produced; only the overview data (real EO, 120% RMT), computed
# with the same lines as in the page loop of the original (load(s, "EO"), trials at 120% RMT, ROI mean +- SEM across trials).
for s, main in order:
    d = load(s, "EO")
    if d is None: continue
    m = d["inten"] == 120; Y = d["roi"][:, m]; mu = Y.mean(1); se = Y.std(1, ddof=1) / np.sqrt(m.sum()); ov[s] = (d["t"], mu, se, int(m.sum()))
ks = sorted(ov); nc = 5; nr = int(np.ceil(len(ks) / nc))
fig, ax = plt.subplots(nr, nc, figsize=(7.0, 1.25 * nr), sharex=True, sharey=True)
for i, a in enumerate(ax.flat):
    if i >= len(ks): a.axis("off"); continue
    s = ks[i]; t, mu, se, n = ov[s]; a.axvspan(-5, 10, color="0.85", lw=0) if H.MODE == "repro" else a.axvspan(-5, 10, color="k", lw=0, zorder=3); a.axhline(0, color="0.6", lw=.4)   # rev1108: black
    a.fill_between(t, mu - se, mu + se, color=COL[120], alpha=.3, lw=0); a.plot(t, mu, color=COL[120], lw=.8)
    a.set_title(f"{s:03d}" + ("" if s in io.MAIN17 else " †"), fontsize=6.5); a.set_xlim(-100, 400)
for a in ax[-1]: a.set_xlabel("Time (ms)")
import os as _o3; ALT = _o3.environ.get("LBL_ALT") == "1"
for a in ax[:, 0]: a.set_ylabel("µV" if not ALT else "Potential (µV)")   # 1108n: alt only
fig.suptitle("Real stimulation, EO, 120% RMT: ROI TEP (0.1–80 Hz), mean ± SEM across trials; † not in the main sample" if CMODE == "repro" else
             ("Real stimulation, EO, 120% RMT: ROI TEP (0.1–80 Hz), mean ± SEM across trials, † not in the main sample" if _o3.environ.get("T_MODE", "rev") == "repro" else
              "Real stimulation, EO, 120% RMT: ROI TEP (0.1–80 Hz), mean ± SE across trials, † not in the main sample"), fontsize=7)   # 1108t: SEM -> SE   # rev1108c: ; -> ,
fig.tight_layout(rect=(0, 0, 1, .97))
# rev1108c: output to Figure_S2_2026-11-08 (nothing written to Figures_all_2026-11-08); same save procedure as rev1108_helpers.save_rev.
# REV1108C_MODE=repro: the rev1108 drawing unchanged -> _check/repro_FigureS8_rev.png (pixel check against Figures_all_2026-11-08/FigureS8_rev.png); rev: commas.
from PIL import Image as _Im
_D = Path(data_dir + r"\Revision_2026\Manuscript_Figures\Figures_2026-11-08t")   # 1108t
assert CMODE == "rev" and ALT; _st = "FigureS3_1108t" if _o3.environ.get("T_MODE", "rev") == "rev" else "repro_FigureS3_1108n"; _d = _D if _o3.environ.get("T_MODE", "rev") == "rev" else _D / "_check" / "_scratch"   # 1108t
fig.savefig(_d / f"{_st}.png", dpi=600, facecolor="white"); fig.savefig(_d / f"{_st}.pdf", facecolor="white")
_im = _Im.open(_d / f"{_st}.png").convert("RGB"); _im.save(_d / f"{_st}.png", dpi=(600, 600))

print("saved", len(ks), "participants in the overview; others:", others)
