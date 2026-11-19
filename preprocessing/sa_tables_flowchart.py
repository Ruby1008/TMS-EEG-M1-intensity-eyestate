"""S-A completion (user 2026-10-18, section 4.3): preprocessing flowchart, ICA statistics real vs sham, and a check of Table S2 against
the final parameters of HANDOFF section 1-1 (Revision_2026\\Preproc\\HANDOFF_20261012.md).
Inputs: Preproc\\QC\\QC_summary_wide.csv (all processed recordings), Preproc\\ica_overrides.csv (manual ICA decisions).
Outputs: data/SA_ICA_stats.csv, data/SA_parameter_check.csv, FigS_preprocessing_flowchart.png/.tiff (600 dpi)."""
import textwrap, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import tep_io as io

P = io.REV / "Preproc"
Q = pd.read_csv(P / "QC" / "QC_summary_wide.csv", encoding="utf-8-sig"); Q = Q[Q.status == "processed"]
O = pd.read_csv(P / "ica_overrides.csv", comment="#", encoding="utf-8-sig"); O = O.dropna(subset=["round"]); O["stim"] = O.cond.str.split("_").str[0]
COLS = [("R1_total", "ICA 1: components removed"), ("R1_decay", "ICA 1: decay rule"), ("R1_tmsMuscle_only", "ICA 1: TMS-muscle classifier only"), ("R1_varPct", "ICA 1: variance removed (%)"),
        ("R2_total", "ICA 2: components removed"), ("R2_blink", "ICA 2: blink"), ("R2_eyeMove", "ICA 2: eye movement"), ("R2_muscle", "ICA 2: persistent muscle"),
        ("R2_elecNoise", "ICA 2: electrode noise"), ("R2_varPct", "ICA 2: variance removed (%)"), ("nManualOverrides", "manual decisions applied")]
rows = []
for col, lab in COLS:
    r = {"measure": lab, "column": col}
    for stim in ("Real", "Sham"):
        for e in ("EO", "EC"):
            v = Q[(Q.stim == stim) & (Q.eyes == e)][col].astype(float)
            r[f"{stim} {e}"] = (f"{np.median(v):.1f}, {v.min():.1f}–{v.max():.1f}" if "Pct" in col else f"{np.median(v):.0f}, {v.min():.0f}–{v.max():.0f}") + f" (n = {len(v)})"
        v = Q[Q.stim == stim][col].astype(float); r[f"{stim} total"] = f"{v.sum():.0f}" if "Pct" not in col else ""
    rows.append(r)
for rd in (1, 2):
    for act in ("remove", "keep"):
        rows.append({"measure": f"manual decision, ICA {rd}, {act}", "column": "ica_overrides.csv", **{f"{stim} total": int(((O.stim == stim) & (O["round"] == rd) & (O.action == act)).sum()) for stim in ("Real", "Sham")}})
S = pd.DataFrame(rows); S.to_csv(io.DATA / "SA_ICA_stats.csv", index=False); print(S.to_string(index=False))
# ---- Table S2 vs HANDOFF 1-1 (final parameters)
chk = [("Pulse detection", "pop_tesa_findpulse FC3, refract 3, rate Real 1000 / Sham 300; repeated detections within 4 s: first kept; epochs locked to the pulse", "yes", ""),
       ("Intensity labels", "trigger code within 100 ms before the pulse; triallist for 6 recordings; verified by pulse size; label_overrides.csv: 028_Real_EC trial 8 -> 120%", "partly", "add the single manual label correction"),
       ("Epoch / first baseline", "-2..2 s; baseline -1000..-10 ms", "yes", ""), ("Removal window", "[-5 10] ms", "yes", ""),
       ("Bad channels", "robust z > 5 of pre-stimulus log-variance; spherical interpolation at the end", "yes", ""), ("Bad trials", "pop_jointprob 5 SD local and global", "yes", ""),
       ("ICA 1", "PCA 40; FastICA symm, tanh, rng 1; decay rule A u C (C: 13-30 ms >= 25x pre [-500 -10] and peak >= 10; A: ratio >= 15, max <= 16 ms, half by 40 ms) + tmsMuscle [11 30] threshold 8", "yes", "rule details in the main text"),
       ("Interpolation and filter", "cubic -> 0.1-80 Hz + 59-61 Hz band-stop (4th-order Butterworth, zero phase) -> cubic", "yes", ""),
       ("ICA 2", "FastICA without PCA; blink Fp1/Fp2 2.5; eyeMove F7/F8 2; muscle slope -0.31 7-70 Hz; elecNoise z 6; tmsMuscle off; fallback electrodes (Fp1/Fp2 bad -> Fpz/AF3/AF4; F7/F8 bad -> F5/F6 or FT7/FT8)", "partly", "add: TMS-muscle classifier off in ICA 2; fallback electrodes"),
       ("Protection rule (ICA 2)", "elecNoise/muscle only: variance > 10% and (ROI back-projection >= 30% of estimated TEP or |z| >= 6 at ROI) -> kept and flagged for manual decision", "yes", ""),
       ("Manual ICA decisions", "ica_overrides.csv, 74 rows (Real 39, Sham 35)", "no", "add row"),
       ("Final steps", "linear interpolation of the removal window -> EOG removed -> interpolate bad channels -> average reference -> baseline -100..-6 ms", "yes", ""),
       ("Outputs", "_F_ERM (no further low-pass), _F_TEP (30 Hz 4th-order Butterworth zero phase)", "yes", ""),
       ("Condition exclusion", "residual 11-25 ms > 50 µV -> exclude condition; 20-50 µV -> sensitivity", "yes", "")]
pd.DataFrame(chk, columns=["step", "final parameter (HANDOFF 1-1)", "in Table S2", "action"]).to_csv(io.DATA / "SA_parameter_check.csv", index=False)
# ---- flowchart
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 6.5})
steps = ["Raw recording (.cnt, 64 channels incl. EOG and mastoids)", "TMS pulse detection on FC3 (threshold: real 1000, sham 300 µV/ms) and intensity labels",
         "Epochs −2 to 2 s locked to the pulse; baseline −1000 to −10 ms", "Pulse window −5 to 10 ms removed", "Bad channels (robust z > 5) and bad epochs (joint probability, 5 SD)",
         "ICA 1 (PCA 40, FastICA): decay rule + TMS-muscle classifier; manual review", "Cubic interpolation; 0.1–80 Hz band-pass, 59–61 Hz band-stop; cubic interpolation",
         "ICA 2 (FastICA): blink, eye movement, persistent muscle, electrode noise; protection rule; manual review",
         "Linear interpolation; EOG removed; channels interpolated; average reference; baseline −100 to −6 ms",
         "Condition exclusion: residual 11–25 ms > 50 µV (20–50 µV: sensitivity analysis)"]
outs = [("_F_ERM (0.1–80 Hz, 61 scalp channels)", "First-layer CEEMDAN (IMFs) → ERM (IMFs 4–7), ITPC, HHCFPC (second layer); single-participant TEPs (S-D); alpha power"),
        ("_F_TEP (30 Hz low-pass)", "Conventional TEP; pre-stimulus statistics (S-B); reliable-response analysis real vs sham (S-F)")]
fig, ax = plt.subplots(figsize=(6.3, 6.4)); ax.set_xlim(0, 10); ax.set_ylim(.6, len(steps) + 3.0); ax.axis("off")
y = len(steps) + 2.6
for k, s in enumerate(steps):
    yy = y - k; ax.add_patch(FancyBboxPatch((1.0, yy - .32), 8.0, .64, boxstyle="round,pad=0.02,rounding_size=0.12", fc="#F2F2F2", ec="0.35", lw=.7))
    ax.text(5.0, yy, textwrap.fill(f"{k + 1}. {s}", 92), ha="center", va="center")
    ax.annotate("", xy=(5.0, yy - .34 - .02), xytext=(5.0, yy - .66 + .02), arrowprops=dict(arrowstyle="<-", lw=.7, color="0.35"))
yb = y - len(steps) - .25
for j, (h, d) in enumerate(outs):
    x0 = .3 + j * 4.85; ax.add_patch(FancyBboxPatch((x0, yb - .75), 4.55, 1.2, boxstyle="round,pad=0.02,rounding_size=0.12", fc="white", ec="0.35", lw=.8))
    ax.text(x0 + 2.275, yb + .2, h, ha="center", va="center", weight="bold"); ax.text(x0 + 2.275, yb - .3, textwrap.fill(d, 58), ha="center", va="center", fontsize=6)
ax.text(5.0, .7, "Real and sham recordings: identical pipeline except the pulse-detection threshold.", ha="center", va="bottom", style="italic")
fig.tight_layout()
for ext in ("png", "tiff"): fig.savefig(io.ROOT / f"FigS_preprocessing_flowchart.{ext}", dpi=600, **({"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tiff" else {}))
print("saved")
