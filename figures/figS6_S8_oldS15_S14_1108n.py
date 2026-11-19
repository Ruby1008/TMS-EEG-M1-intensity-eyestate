# _1108n copy (user 2026-10-28, Figures_2026-11-08n) of Figures_all_2026-11-08/code/figS14_S15_r21_n100_confound_rev1108.py. Changes: N100 dashed lines and label removed (B); Cz -> FigureS6_1108n, ROI -> FigureS8_1108n; LBL_ALT=1: y label 'Potential (µV)' -> _check/alt_potential.
# _rev1108 copy (user 2026-10-27, Figures_all_2026-11-08) of Manuscript_Figures/code/r21_n100_confound.py (FigS_R21_ERM_real_sham_ROI.png = Figure S14,
# FigS_R21_ERM_real_sham_Cz.png = Figure S15, MD5). Statistics recomputed in memory only (no CSV written; the sliding-window q values are asserted equal to
# data/R21_diffwave_slidingwindow.csv). Change (REV1108_MODE=rev): grey band 120-140 ms -> two black 0.8-pt dashed lines with the label 'N100' at the top of
# each panel (S14: just below the blue significance bars at 4.5 µV; S15: right-aligned at 140 ms). Outputs FigureS14_rev / FigureS15_rev (.png 600 dpi, .pdf).
"""N100 confound (master plan I): real vs sham ERM (CEEMDAN IMF4-7, first layer, same parameters) at the ROI and at Cz, and a sliding-window 2 (eye) x 3 (intensity)
rmANOVA on the real - sham difference wave (21-sample windows, 10-ms steps, 0-400 ms, Greenhouse-Geisser, BH-FDR across windows per effect).
Participants: the 17 main participants with sham data in both eye states (n = 16; 024 sham EO missing). ERM files: HHSA\\erm_spTMS_<sub>_<Real|Sham>_<eye>.mat.
Also the N100 (120-140 ms) real vs sham per condition (paired t, Holm across six).
Outputs: data/R21_diffwave_slidingwindow.csv, data/R21_N100_real_vs_sham.csv, FigS_R21_ERM_real_sham_ROI.png/.tiff, FigS_R21_ERM_real_sham_Cz.png/.tiff (600 dpi)
"""
data_dir = r"path/to/data"   # release: set this path
import numpy as np, pandas as pd, pingouin as pg, scipy.io as sio
from scipy import stats
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
    new = {'FigureS14_rev': 'FigureS8_1108n', 'FigureS15_rev': 'FigureS6_1108n'}[stem]
    if new is None: return
    if ALT: H.OUT = N_OUT / "_check" / "alt_potential"; new = new + "_alt"
    _save0(fig, new, **kw); H.OUT = N_OUT; print("saved", new)
H.save_rev = _save1108n
WL = []

def erm(s, stim, e): return sio.loadmat(io.HHSA / f"erm_spTMS_{s:03d}_{stim}_{e}.mat", squeeze_me=True)
SUBS = [s for s in io.MAIN17 if all((io.HHSA / f"erm_spTMS_{s:03d}_Sham_{e}.mat").exists() for e in io.EYES)]
R0 = erm(SUBS[0], "Real", "EO"); t = np.round(R0["time"]).astype(int); labs = [str(x) for x in R0["chanlabels"]]; cz = labs.index("Cz")
W = {}
for s in SUBS:
    for stim in ("Real", "Sham"):
        for e in io.EYES:
            R = erm(s, stim, e); W[(s, stim, e, "ROI")] = R["erm_roi"]; W[(s, stim, e, "Cz")] = R["erm_allch"][cz]       # 3 x time
starts = np.arange(0, 400, 10)
def bh(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); q = np.empty(m); q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]; return np.minimum(q, 1)
rows = []
for a0 in starts:
    i0 = np.where(t >= a0)[0][0]
    df = pd.DataFrame([dict(sub=s, eye=e, intensity=lv, y=(W[(s, "Real", e, "ROI")][k] - W[(s, "Sham", e, "ROI")][k])[i0:i0 + 21].mean()) for s in SUBS for e in io.EYES for k, lv in enumerate(io.LEV)])
    a = pg.rm_anova(data=df, dv="y", within=["eye", "intensity"], subject="sub", correction=True, effsize="np2").set_index("Source"); r = dict(win_start_ms=a0)
    for src, k in (("eye", "eye"), ("intensity", "intensity"), ("eye * intensity", "interaction")):
        x = a.loc[src]; r.update({f"F_{k}": x["F"], f"p_{k}": x["p_GG_corr"] if pd.notna(x.get("p_GG_corr", np.nan)) else x["p_unc"], f"np2_{k}": x["np2"]})
    rows.append(r)
R = pd.DataFrame(rows)
for k in ("eye", "intensity", "interaction"): R[f"q_{k}"] = bh(R[f"p_{k}"])
_R0 = pd.read_csv(io.DATA / "R21_diffwave_slidingwindow.csv"); assert np.allclose(_R0.q_intensity.values, R.q_intensity.values)   # rev1108: not written; checked against the stored file
w = (t >= 120) & (t <= 140); N = []
for e in io.EYES:
    for k, lv in enumerate(io.LEV):
        a = np.array([W[(s, "Real", e, "ROI")][k][w].mean() for s in SUBS]); b = np.array([W[(s, "Sham", e, "ROI")][k][w].mean() for s in SUBS]); r = stats.ttest_rel(a, b)
        N.append(dict(eye=e, intensity=lv, n=len(SUBS), real=a.mean(), sham=b.mean(), t=r.statistic, p=r.pvalue))
N = pd.DataFrame(N); N["p_holm"] = pg.multicomp(N.p.values, method="holm")[1]   # rev1108: data/R21_N100_real_vs_sham.csv not written
def runs(m):
    out, st, prev = [], None, None
    for a0, x in zip(starts, m):
        if x and st is None: st = a0
        if not x and st is not None: out.append(f"{st}-{prev + 20}"); st = None
        prev = a0
    if st is not None: out.append(f"{st}-{prev + 20}")
    return ", ".join(out) or "none"
print("n =", len(SUBS), "| diff-wave intensity FDR:", runs(R.q_intensity < .05), "| eye FDR:", runs(R.q_eye < .05), "| interaction FDR:", runs(R.q_interaction < .05))
pd.set_option("display.width", 200); print(N.round(4).to_string(index=False))
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 9})
for site in ("ROI", "Cz"):
    fig, axs = plt.subplots(2, 3, figsize=(11, 5.6), sharex=True, sharey=True)
    for r_, e in enumerate(io.EYES):
        for c_, lv in enumerate(io.LEV):
            ax = axs[r_, c_]; k = io.LEV.index(lv)
            for stim, col in (("Real", "#c0392b"), ("Sham", "0.3")): ax.plot(t, np.mean([W[(s, stim, e, site)][k] for s in SUBS], 0), color=col, lw=1.3, label=stim)
            ax.plot(t, np.mean([W[(s, "Real", e, site)][k] - W[(s, "Sham", e, site)][k] for s in SUBS], 0), color="#2a78d6", lw=1, ls="--", label="Real − Sham")
            if site == "ROI":
                for a0 in starts[(R.q_intensity < .05).values]: ax.plot([a0, a0 + 20], [ax.get_ylim()[1] if False else 4.5] * 2, color="#2a78d6", lw=3)
            if H.MODE == "repro": ax.axvspan(120, 140, color="0.92", lw=0, zorder=0)
            elif False: WL.append((site, ax, H.dashed_window(ax, 120, 140, "N100", **(dict(y_data=4.5, pad_pt=3.5) if site == "ROI" else dict(ha="right", x=140)))))   # rev1108: dashed lines + label "N100" at the top
            # of the panel (ROI: just below the row of blue bars at 4.5 µV, which the label would otherwise cover; Cz: right-aligned at 140 ms, centred it
            # touches the sham curve of EC 120% RMT)
            ax.axhline(0, color="0.6", lw=.5); ax.set_xlim(0, 400); ax.set_title(f"{e}, {lv}% RMT")
            if c_ == 0: ax.set_ylabel(f"ERM, {site} (µV)" if not ALT else "Potential (µV)")   # 1108n: alt version only (DECISION NEEDED)
            if r_ == 1: ax.set_xlabel("Time (ms)")
    axs[0, 0].legend(frameon=False, fontsize=7)
    fig.suptitle(f"ERM (IMF4–7) real vs sham stimulation, {site}", fontsize=9)
    fig.tight_layout(); k_ = {"ROI": "FigS14", "Cz": "FigS15"}[site]; H.save_signature(fig, k_)
    if False:   # label clearance (no window labels in 1108n): curves and the blue 2-point bars (data coordinates)
        import json as _j
        cl = [dict(panel=a_.get_title(), **H.label_clearance(fig, a_, t_, [(l.get_xdata(), l.get_ydata()) for l in a_.lines if l.get_gid() != "rev1108_dash" and l.get_transform() == a_.transData and np.size(l.get_xdata()) >= 2])) for s_, a_, t_ in WL if s_ == site]
        _j.dump(cl, open(H.CK / f"{k_.replace('FigS', 'FigureS')}_window_labels.json", "w"), indent=1)
    if H.MODE == "repro": fig.savefig(H.CK / f"repro_{k_}.png", dpi=600)
    else: H.save_rev(fig, k_.replace("FigS", "FigureS") + "_rev")
    plt.close(fig)
