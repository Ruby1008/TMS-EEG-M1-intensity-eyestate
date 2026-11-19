# _1108n copy (user 2026-10-28, Figures_2026-11-08n) of Figures_all_2026-11-08/code/figS6_S7_a3_rev1108.py. Changes: none in the main output (A: DECISION NEEDED); LBL_ALT=1: old S7 (new S4) y label 'Potential (µV)' -> _check/alt_potential/FigureS4_1108n_alt.png; old S6 not saved.
# _rev1108 copy (user 2026-10-27, Figures_all_2026-11-08) of qc/figures_a3/figS6_S7_a3.py (FigS6_a3.png = Figure S6, FigS7_a3.png = Figure S7, MD5).
# Changes (REV1108_MODE=rev): S6 dark-grey removal bar (black, alpha .8) -> black; S7 dark-grey removal bar (alpha .75) -> black, grey band -20..-5 ms ->
# two black 0.8-pt dashed lines with the label '−20 to −5 ms' at the top of each panel (right-aligned at -6 ms, left of the black bar), title 'gray band' ->
# 'dashed'. Data and everything else unchanged. Outputs FigureS6_rev / FigureS7_rev (.png 600 dpi, .pdf).
"""Part 4 A3 (user 2026-10-22): Supplementary Figures S6 and S7 redrawn in the layout of the originals.
Sources (found in part 4; round-15 audit had listed them as missing): MF/code/erm_vs_tep_prestim.py (SUBSET=main17 -> FigS_ERMvsTEP_overlay_main17.png = Figure S6)
and MF/code/erm_prestim_variants.py (SUBSET=main17 -> FigS_ERMvsTEP_variants_main17.png = Figure S7). Only the plotting parts are reproduced here, with the
same data and sample definitions (HHSA/erm_*.mat, data/TEP_roi_waveforms.csv, Preproc/QC/QC_summary_long.csv stage D, MAIN17); no analysis output is written.
Changes: Times New Roman, line widths as the other line figures (1.5 main trace, 1.2 others), every line defined in the legend (the S6 legend had dropped the
two black lines because their labels began with '_'), en dashes and 'gray'. Outputs: qc/figures_a3/FigS6_a3.png/.tiff, FigS7_a3.png/.tiff (600 dpi)."""
data_dir = r"path/to/data"   # release: set this path
import sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.io as sio
from scipy import signal
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
MF = Path(data_dir + r"\Revision_2026\Manuscript_Figures"); sys.path.insert(0, str(MF / "code")); from tep_io import MAIN17
sys.path.insert(0, data_dir + r"\Revision_2026\Manuscript_Figures\Figures_all_2026-11-08\code"); import rev1108_helpers as RH

# ---- 1108n (user 2026-10-28): outputs only to Figures_2026-11-08n; the shared helper is redirected so nothing is written to Figures_all_2026-11-08.
import os as _os1, json as _js1
from pathlib import Path as _P1
N_OUT = _P1(data_dir + r"\Revision_2026\Manuscript_Figures\Figures_2026-11-08n"); ALT = _os1.environ.get("LBL_ALT") == "1"
RH.OUT = N_OUT; RH.CK = N_OUT / "_check" / "_scratch"; RH.CK.mkdir(parents=True, exist_ok=True); RH.save_signature = lambda *a, **k: None
_save0 = RH.save_rev
def _save1108n(fig, stem, **kw):
    new = {'FigureS6_rev': None, 'FigureS7_rev': 'FigureS4_1108n'}[stem]
    if new is None: return
    if ALT: RH.OUT = N_OUT / "_check" / "alt_potential"; new = new + "_alt"
    _save0(fig, new, **kw); RH.OUT = N_OUT; print("saved", new)
RH.save_rev = _save1108n
OUT = Path(__file__).parent; HH = MF / "HHSA"
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 10, "axes.linewidth": 0.8, "mathtext.fontset": "stix"})
Q = pd.read_csv(MF.parent / "Preproc" / "QC" / "QC_summary_long.csv", encoding="utf-8-sig")[["id", "intensity", "stageD"]]
W = pd.read_csv(MF / "data" / "TEP_roi_waveforms.csv"); tcols = [c for c in W.columns if c.startswith("t_")]; tt = np.arange(-2000, -2000 + len(tcols))
LEV = [80, 100, 120]
def save(fig, name):   # rev1108: repro -> _check/repro_<FigS6|FigS7>.png (original call), rev -> Figures_all_2026-11-08/FigureS6_rev / FigureS7_rev (.png 600 dpi, .pdf)
    k = {"FigS6_a3": "FigS6", "FigS7_a3": "FigS7"}[name]; RH.save_signature(fig, k)
    if RH.MODE == "repro": fig.savefig(RH.CK / f"repro_{k}.png", dpi=600)
    else: RH.save_rev(fig, k.replace("FigS", "FigureS") + "_rev")
    plt.close(fig)
# ---------------- Figure S6 (erm_vs_tep_prestim.py, panel (a)); sample: id with any non-excluded intensity (stage D) among the 17 main participants
wave = {}
for f in sorted(HH.glob("erm_spTMS_*.mat")):
    R = sio.loadmat(f, squeeze_me=True); wave[str(R["id"])] = (R["time"], R["erm_roi"])
main_ids = {i for i in Q[Q.stageD != "exclude"].id.unique() if int(i[6:9]) in MAIN17}
fig, axs = plt.subplots(2, 3, figsize=(13, 6), sharex=True, sharey=True)
cols = {80: "#f9954e", 100: "#f97213", 120: "#c55505"}
for r_, e in enumerate(("EO", "EC")):
    for c_, lv in enumerate(LEV):
        ax = axs[r_, c_]; ids_ = [i for i in wave if i.endswith(f"Real_{e}") and i in main_ids]
        k = LEV.index(lv); t = wave[ids_[0]][0]; erm = np.mean([wave[i][1][k] for i in ids_], axis=0)
        ax.axvspan(-5, 10, color="k", alpha=.8, lw=0) if RH.MODE == "repro" else ax.axvspan(-5, 10, color="k", lw=0, zorder=3); ax.axhline(0, color="#999", lw=.7)   # rev1108: black
        ax.plot(t, erm, color=cols[lv], lw=1.5, label="ERM (CEEMDAN IMF4–7)")
        for fl, ls, lab in (("TEP", "-", "Conventional TEP (0.1–30 Hz)"), ("ERM", ":", "0.1–80 Hz data without the 30 Hz low-pass")):
            g = W[(W.file == fl) & (W.intensity == lv) & W.id.isin(ids_)]
            ax.plot(tt, g[tcols].to_numpy().mean(0), color="k", lw=1.2, ls=ls, label=lab)
        ax.set_xlim(-100, 400); ax.set_title(f"Real {e} {lv}% RMT", fontsize=11)
        if r_ == 0 and c_ == 0: ax.legend(fontsize=8, loc="lower right", frameon=False)
for ax in axs[1]: ax.set_xlabel("Time (ms)")
for ax in axs[:, 0]: ax.set_ylabel("ROI (µV)")
fig.suptitle("Grand average, ROI FC3/FC1/C3/C1: ERM (CEEMDAN IMF4–7) vs conventional TEP", fontsize=12); fig.tight_layout(); save(fig, "FigS6_a3")
# ---------------- Figure S7 (erm_prestim_variants.py); sample: id x intensity not excluded at stage D, 17 main participants
b, a = signal.butter(4, [4.2, 36.4], btype="band", fs=1000); waves = {}
for f in sorted(HH.glob("erm_spTMS_*_Real_*.mat")):
    R = sio.loadmat(f, squeeze_me=True); id_ = str(R["id"]); t = R["time"]; erm = R["erm_roi"]
    for k, lv in enumerate(LEV):
        g = W[(W.id == id_) & (W.intensity == lv)]
        if not len(g): continue
        x_erm = g[g.file == "ERM"][tcols].to_numpy().ravel(); x_tep = g[g.file == "TEP"][tcols].to_numpy().ravel()
        if Q[(Q.id == id_) & (Q.intensity == lv)].stageD.iloc[0] != "exclude" and int(id_[6:9]) in MAIN17:
            w = erm[k]; waves.setdefault((id_[-2:], lv), []).append((w - w[(t >= -300) & (t <= -100)].mean(), signal.filtfilt(b, a, x_erm), x_tep))
t = sio.loadmat(next(HH.glob("erm_spTMS_*_Real_*.mat")), squeeze_me=True)["time"]
fig, axs = plt.subplots(2, 3, figsize=(14, 6.5), sharex=True); WL7 = []
for r_, e in enumerate(("EO", "EC")):
    for c_, lv in enumerate(LEV):
        ax = axs[r_, c_]; L = waves[(e, lv)]
        erm = np.mean([x[0] for x in L], 0); bp = np.mean([x[1] for x in L], 0); tp = np.mean([x[2] for x in L], 0)
        if RH.MODE == "repro": ax.axvspan(-5, 10, color="k", alpha=.75, lw=0); ax.axvspan(-20, -5, color="#ddd", lw=0)
        else: ax.axvspan(-5, 10, color="k", lw=0, zorder=3); WL7.append((ax, RH.dashed_window(ax, -20, -5, "−20 to −5 ms", fontsize=9, ha="right", x=-6)))   # rev1108
        ax.axhline(0, color="#999", lw=.7)
        h1, = ax.plot(t, erm, color="#c55505", lw=1.5); h2, = ax.plot(t, erm - erm[(t >= -100) & (t <= -10)].mean(), color="#f9954e", lw=1.2, ls="--")
        h3, = ax.plot(tt, bp - bp[(tt >= -300) & (tt <= -100)].mean(), color="#2a78d6", lw=1.2); h4, = ax.plot(tt, tp - tp[(tt >= -300) & (tt <= -100)].mean(), color="k", lw=1.2)
        H = [h1, h2, h3, h4]
        ax.set_xlim(-150, 300); ax.set_title(f"Real {e} {lv}% RMT", fontsize=11)
for ax in axs[1]: ax.set_xlabel("Time (ms)")
for ax in axs[:, 0]: ax.set_ylabel("ROI (µV)" if not ALT else "Potential (µV)")   # 1108n: alt only
fig.legend(H, ["ERM IMF4–7, reference −300 to −100 ms", "ERM, baseline −100 to −10 ms", "0.1–80 Hz data band-passed 4.2–36.4 Hz (no EMD)", "Conventional TEP (0.1–30 Hz)"], loc="lower center", ncol=4, fontsize=10, frameon=False)
fig.suptitle("Grand averages (ROI), gray band: pre-stimulus window −20 to −5 ms" if RH.MODE == "repro" else "Grand averages (ROI), dashed: pre-stimulus window −20 to −5 ms", fontsize=12); fig.tight_layout(rect=(0, 0.06, 1, 1))
if False:   # label clearance (curves of each panel)
    import json as _j; _j.dump([dict(panel=a_.get_title(), **RH.label_clearance(fig, a_, t_, [(l.get_xdata(), l.get_ydata()) for l in a_.lines if np.size(l.get_xdata()) > 2 and l.get_gid() != "rev1108_dash"])) for a_, t_ in WL7], open(RH.CK / "FigureS7_window_labels.json", "w"), indent=1)
save(fig, "FigS7_a3")
print({k: len(v) for k, v in waves.items()}, {e: len([i for i in wave if i.endswith(f"Real_{e}") and i in main_ids]) for e in ("EO", "EC")})
