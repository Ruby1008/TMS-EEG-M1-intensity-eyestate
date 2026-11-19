# _1108n copy (user 2026-10-28, Figures_2026-11-08n) of Figures_all_2026-11-08/code/figS21_aep_compare_erm_itpc_rev1108.py. Changes: N100 dashed lines and label removed (B); output FigureS12_1108n; LBL_ALT=1: y label 'Potential (µV)'.
# _rev1108 copy (user 2026-10-27, Figures_all_2026-11-08) of Manuscript_Figures/code/aep_compare_erm_itpc.py (FigS3_ERM_AEP_before_after.png = Figure S21, MD5).
# ERM statistics recomputed in memory only (no CSV written; q values asserted equal to data/AEP_ERM_slidingwindow_before_after.csv); ITPC overlap part removed.
# Change (REV1108_MODE=rev): grey band 120-140 ms -> two black 0.8-pt dashed lines with the label 'N100' at the top of each panel, just below the two rows
# of significance bars. Output FigureS21_rev.png/.pdf.
"""AEP-ICA robustness (master plan J, S3-S5): ERM and ITPC before vs after AEP removal, 17 participants.
Before = HHSA\\erm_*.mat (main analysis); after = HHSA_AEPr\\erm_*.mat (aep_remove.m -> run_layer1_aepr.m, same CEEMDAN parameters/seed).
ERM: sliding-window 2 x 3 rmANOVA (21 samples, 10-ms steps, 0-400 ms, GG, BH-FDR per effect) and N100 (120-140 ms) ANOVA, before and after.
ITPC: Fig 7 analysis (itpc_fig78_G.py; after = ITPC_VARIANT=AEPr) -> overlap of the FDR masks per effect and IMF.
P30/N45/P60: after the component windows are decided.
Outputs: data/AEP_ERM_slidingwindow_before_after.csv, data/AEP_ERM_N100_anova.csv, data/AEP_ITPC_mask_overlap.csv,
         FigS3_ERM_AEP_before_after.png/.tiff (600 dpi)
"""
data_dir = r"path/to/data"   # release: set this path
import numpy as np, pandas as pd, pingouin as pg
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import sys as _s0; _s0.path.insert(0, data_dir + r'\Revision_2026\Manuscript_Figures\code')
import tep_io as io
import sys as _s; _s.path.insert(0, data_dir + r'\Revision_2026\Manuscript_Figures\Figures_all_2026-11-08\code'); import rev1108_helpers as H

# ---- 1108n (user 2026-10-28): outputs only to Figures_2026-11-08n; the shared helper is redirected so nothing is written to Figures_all_2026-11-08.
import os as _os1, json as _js1
from pathlib import Path as _P1
N_OUT = _P1(data_dir + r"\Revision_2026\Manuscript_Figures\Figures_2026-11-08n"); ALT = _os1.environ.get("LBL_ALT") == "1"
H.OUT = N_OUT; H.CK = N_OUT / "_check" / "_scratch"; H.CK.mkdir(parents=True, exist_ok=True); H.save_signature = lambda *a, **k: None
_save0 = H.save_rev
def _save1108n(fig, stem, **kw):
    new = {'FigureS21_rev': 'FigureS12_1108n'}[stem]
    if new is None: return
    if ALT: H.OUT = N_OUT / "_check" / "alt_potential"; new = new + "_alt"
    _save0(fig, new, **kw); H.OUT = N_OUT; print("saved", new)
H.save_rev = _save1108n
WL = []

S17, LEV, EYES = io.MAIN17, io.LEV, io.EYES
FOLD = {"before": io.HHSA, "after": io.ROOT / "HHSA_AEPr"}
ERM = {v: {(s, e): io.load_erm(s, e, f) for s in S17 for e in EYES} for v, f in FOLD.items()}
t = np.round(ERM["before"][(S17[0], "EO")]["time"]).astype(int); starts = np.arange(0, 400, 10)
def bh(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); q = np.empty(m); q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]; return np.minimum(q, 1)
def anova(df):
    a = pg.rm_anova(data=df, dv="y", within=["eye", "intensity"], subject="sub", correction=True, effsize="np2").set_index("Source"); out = {}
    for src, k in (("eye", "eye"), ("intensity", "intensity"), ("eye * intensity", "interaction")):
        r = a.loc[src]; out.update({f"F_{k}": r["F"], f"p_{k}": r["p_GG_corr"] if pd.notna(r.get("p_GG_corr", np.nan)) else r["p_unc"], f"np2_{k}": r["np2"],
                                    f"eps_{k}": r["eps"] if pd.notna(r.get("eps", np.nan)) else 1.0})
    return out
SW, N1 = [], []
for v in FOLD:
    rows = []
    for a0 in starts:
        i0 = np.where(t >= a0)[0][0]
        df = pd.DataFrame([dict(sub=s, eye=e, intensity=lv, y=ERM[v][(s, e)]["erm_roi"][k][i0:i0 + 21].mean()) for s in S17 for e in EYES for k, lv in enumerate(LEV)])
        rows.append(dict(version=v, win_start_ms=a0, **anova(df)))
    R = pd.DataFrame(rows)
    for k in ("eye", "intensity", "interaction"): R[f"q_{k}"] = bh(R[f"p_{k}"])
    SW.append(R)
    w = (t >= 120) & (t <= 140)
    df = pd.DataFrame([dict(sub=s, eye=e, intensity=lv, y=ERM[v][(s, e)]["erm_roi"][k][w].mean()) for s in S17 for e in EYES for k, lv in enumerate(LEV)])
    N1.append(dict(version=v, component="N100", window_ms="120-140", n=len(S17), **anova(df)))
SW = pd.concat(SW); _S0 = pd.read_csv(io.DATA / "AEP_ERM_slidingwindow_before_after.csv"); assert np.allclose(_S0.q_intensity.values, SW.q_intensity.values)   # rev1108: not written; checked against the stored file
N1 = pd.DataFrame(N1)   # rev1108: not written
def runs(m):
    out, st = [], None
    for a0, x in zip(starts, m):
        if x and st is None: st = a0
        if not x and st is not None: out.append(f"{st}-{prev + 20}"); st = None
        prev = a0
    if st is not None: out.append(f"{st}-{prev + 20}")
    return ", ".join(out) or "none"
for v in FOLD:
    r = SW[SW.version == v]; print(v, "ERM intensity FDR:", runs(r.q_intensity < .05), "| eye FDR:", runs(r.q_eye < .05), "| interaction FDR:", runs(r.q_interaction < .05))
pd.set_option("display.width", 220); print(N1.round(4).to_string())
# rev1108: the ITPC mask overlap (data/AEP_ITPC_mask_overlap.csv) is not recomputed or written.
# Fig S3
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 9})
fig, axs = plt.subplots(2, 3, figsize=(11, 5.8), sharex=True, sharey=True)
for r_, e in enumerate(EYES):
    for c_, lv in enumerate(LEV):
        ax = axs[r_, c_]; k = LEV.index(lv)
        for v, sty in (("before", dict(color="k", lw=1.3)), ("after", dict(color="#c0392b", lw=1.3, ls="--"))):
            ax.plot(t, np.mean([ERM[v][(s, e)]["erm_roi"][k] for s in S17], 0), label=f"{'Before' if v == 'before' else 'After'} AEP removal", **sty)
        for v, y, col in (("before", 5.5, "k"), ("after", 5.0, "#c0392b")):
            for a0 in starts[(SW[SW.version == v].q_intensity < .05).values]: ax.plot([a0, a0 + 20], [y, y], color=col, lw=3, solid_capstyle="butt")
        if H.MODE == "repro": ax.axvspan(120, 140, color="0.92", lw=0)
        elif False: WL.append((ax, H.dashed_window(ax, 120, 140, "N100", y_data=5.0, pad_pt=3.5)))   # 1108n: removed (B)   # rev1108: dashed lines + "N100" just below the two rows of bars (5.5, 5.0 µV)
        ax.axhline(0, color="0.6", lw=.5); ax.set_xlim(0, 400); ax.set_title(f"Real {e}, {lv}% RMT")
        if c_ == 0: ax.set_ylabel("ERM, ROI (µV)" if not ALT else "Potential (µV)")   # 1108n: alt version only (DECISION NEEDED)
        if r_ == 1: ax.set_xlabel("Time (ms)")
axs[0, 0].legend(frameon=False, fontsize=7, loc="lower right")
fig.suptitle("ERM (IMF4–7) before and after AEP-ICA removal", fontsize=9)
fig.tight_layout(); H.save_signature(fig, "FigS21")
if False:
    import json as _j
    _j.dump([dict(panel=a_.get_title(), **H.label_clearance(fig, a_, t_, [(l.get_xdata(), l.get_ydata()) for l in a_.lines if l.get_gid() != "rev1108_dash" and l.get_transform() == a_.transData and np.size(l.get_xdata()) >= 2])) for a_, t_ in WL], open(H.CK / "FigureS21_window_labels.json", "w"), indent=1)
if H.MODE == "repro": fig.savefig(H.CK / "repro_FigS21.png", dpi=600)
else: H.save_rev(fig, "FigureS21_rev")
print("saved", H.MODE)
