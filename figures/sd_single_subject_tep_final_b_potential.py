# final_b_potential copy (user 2026-10-28) of final_b/code/sd_single_subject_tep_rev1108d_final_b.py: y label 'Amplitude (µV)' -> 'Potential (µV)'; output final_b_potential/.
# _final_b copy (user 2026-10-27) of FigS_SD_2026-11-08/code/sd_single_subject_tep_rev1108d.py run with REV1108D_STYLE=backup (the (b) file): the only change is
# that the -5..10 ms black bar is drawn only in panels with data (same condition as the 'not available' text: d is None); output final_b/.
# _rev1108d copy (user 2026-10-27, FigS_SD_2026-11-08) of Manuscript_Figures/code/sd_single_subject_tep.py: per-participant PDF only. Changes (rev):
# page titles '—' / ';' -> ',', 'grey:' -> 'black bar:'; removal window -5..10 ms grey (0.85) -> black solid bar, zorder 3, in all 8 panels (as Figure S8).
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
import os
from pathlib import Path
# rev1108d: REV1108D_MODE = repro (the original drawing) | rev (default: titles with commas and "black bar:", black removal window as Figure S8).
# REV1108D_STYLE = current (default: this script, as the PDF in Manuscript_Figures, Times New Roman) | backup (the three page-drawing differences of
# code/_backup_20261019/sd_single_subject_tep.py: Arial, "(n = k trials)" in the ROI panel titles, "(not in the 17-participant main sample)"), written to
# _check/alt_backup_style/ for the user's decision. Outputs only in FigS_SD_2026-11-08; the overview figure (Figure S8) is not produced.
MODE = os.environ.get("REV1108D_MODE", "rev"); STYLE = os.environ.get("REV1108D_STYLE", "backup")   # final_b: (b) style by default; assert MODE in ("rev", "repro") and STYLE in ("current", "backup")
OUT = Path(data_dir + r"\Revision_2026\Manuscript_Figures\FigS_SD_2026-11-08"); CK = OUT / "_check"; D_ = CK / "alt_backup_style" if STYLE == "backup" else CK
PDF_ = (OUT if STYLE == "current" else D_) / "FigS_SD_single_participant_TEPs.pdf" if MODE == "rev" else D_ / "repro_FigS_SD_single_participant_TEPs.pdf"
D_.mkdir(parents=True, exist_ok=True)
# final_b: output to FigS_SD_2026-11-08/final_b (the (a) and (b) files are not touched); repro / current style keep the rev1108d paths unused here
FB = OUT / "final_b_potential"; assert MODE == "rev" and STYLE == "backup", "final_b copy: rev mode, backup style only"   # final_b_potential: output folder
D_ = FB / "_check"; D_.mkdir(parents=True, exist_ok=True); PDF_ = FB / "FigS_SD_single_participant_TEPs.pdf"

plt.rcParams.update({"font.family": "Times New Roman" if STYLE == "current" else "Arial", "font.size": 7}); COL = {80: "#9ecae1", 100: "#4292c6", 120: "#08306b"}
others = sorted((set(io.AVAIL["EO"]) | set(io.AVAIL["EC"])) - set(io.MAIN17)); order = [(s, True) for s in io.MAIN17] + [(s, False) for s in others]
def load(s, e):
    if s not in io.AVAIL[e]: return None
    x, labs, t, inten, _ = io.load_set(io.FINAL, f"spTMS_{s:03d}_Real_{e}_F_ERM"); t = np.round(t).astype(int)
    sc = [i for i, l in enumerate(labs) if l not in ("M1", "M2")]; ri = [labs.index(c) for c in io.ROI]; w = (t >= -100) & (t <= 400)
    return dict(t=t[w], roi=x[ri][:, w].astype(float).mean(0), all=x[sc][:, w].astype(float), inten=inten)
rows, ov = [], {}
with PdfPages(PDF_) as pdf:   # rev1108d: output path
    for s, main in order:
        fig, ax = plt.subplots(2, 4, figsize=(10, 4.6), sharex=True)
        for r, e in enumerate(io.EYES):
            d = load(s, e)
            for c, lv in enumerate(io.LEV):
                a = ax[r, c]
                if d is not None: a.axvspan(-5, 10, color="0.85", lw=0) if MODE == "repro" else a.axvspan(-5, 10, color="k", lw=0, zorder=3)   # final_b: no bar in empty panels
                a.axhline(0, color="0.6", lw=.4)
                if d is None: a.text(.5, .5, "not available", transform=a.transAxes, ha="center", color="0.5"); continue
                m = d["inten"] == lv; Y = d["roi"][:, m]; mu = Y.mean(1); se = Y.std(1, ddof=1) / np.sqrt(m.sum())
                a.fill_between(d["t"], mu - se, mu + se, color=COL[lv], alpha=.35, lw=0); a.plot(d["t"], mu, color=COL[lv], lw=.9)
                a.set_title(f"{e}, {lv}% RMT" if STYLE == "current" else f"{e}, {lv}% RMT (n = {int(m.sum())} trials)"); rows.append(dict(sub=s, main17=main, eye=e, intensity=lv, n_trials=int(m.sum())))
                if e == "EO" and lv == 120: ov[s] = (d["t"], mu, se, int(m.sum()))
            a = ax[r, 3]
            if d is not None: a.axvspan(-5, 10, color="0.85", lw=0) if MODE == "repro" else a.axvspan(-5, 10, color="k", lw=0, zorder=3)   # final_b: no bar in empty panels
            if d is not None:
                m = d["inten"] == 120; B = d["all"][:, :, m].mean(2)
                a.plot(d["t"], B.T, color="0.65", lw=.3); a.plot(d["t"], d["roi"][:, m].mean(1), color="k", lw=.9); a.set_title(f"{e}, 120% RMT, 61 channels")
            for a in ax[r]: a.set_xlim(-100, 400)
            ax[r, 0].set_ylabel("Potential (µV)")   # final_b_potential: EEG potential axis label
        for a in ax[1]: a.set_xlabel("Time (ms)")
        _ns = " (not in the main sample)" if STYLE == "current" else " (not in the 17-participant main sample)"
        if MODE == "repro": _tt = f"Participant {s:03d}" + ("" if main else _ns) + " — real stimulation, ROI FC3/FC1/C3/C1, 0.1–80 Hz (before the 30 Hz low-pass); mean ± SEM across trials; grey: −5 to 10 ms removed"
        else: _tt = f"Participant {s:03d}" + ("" if main else _ns) + ", real stimulation, ROI FC3/FC1/C3/C1, 0.1–80 Hz (before the 30 Hz low-pass), mean ± SEM across trials, black bar: −5 to 10 ms removed"   # rev1108d
        fig.suptitle(_tt, fontsize=8)
        fig.tight_layout(rect=(0, 0, 1, .95)); pdf.savefig(fig); plt.close(fig); print(s)
pd.DataFrame(rows).to_csv(D_ / f"ntrials_{MODE}.csv", index=False)   # rev1108d: not written to data/
print("saved", MODE, STYLE, PDF_)
# rev1108d: the overview figure (FigS_SD_overview_EO120 = Figure S8) is not produced (done in Figure_S2_2026-11-08).
