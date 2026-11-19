# _1109i copy (user 2026-10-28, Figures_2026-11-09i) of Figures_2026-11-08n/code/figS2_1108n.py: section labels in the two output boxes renumbered
# (S-A -> S9, S-B -> S10, S-D -> S12); nothing else changed. S2I_MODE=start draws the 1108n figure (start check against FigureS2_1108n.png, into
# _check); S2I_MODE=rev (default) writes FigureS2_1109i.png/.pdf. Text is re-wrapped by the same textwrap.fill(…, 58); box sizes unchanged.
# _1108n copy (user 2026-10-28, Figures_2026-11-08n) of Figure_S2_2026-11-08/code/figS2_sa_tables_flowchart_rev1108c.py: step 7 one line, output box titles and texts (C).
# REV1108B_MODE=repro still draws the rev1108 drawing (start check of the chain); the rev1108c drawing is the comparison base (FigureS2_noHHCFPC_rev_v2.png).
# _rev1108c copy (user 2026-10-27, Figure_S2_2026-11-08) of figS2_sa_tables_flowchart_rev1108b.py: options (b) of the three [DECISION NEEDED] lines (two-line
# splits); outputs FigureS2_noHHCFPC_rev_v2.png/.pdf. repro (REV1108B_MODE=repro) = the rev1108 drawing, as in rev1108b.
# _rev1108b copy (user 2026-10-27, Figure_S2_2026-11-08) of Figures_all_2026-11-08/code/figS2_sa_tables_flowchart_rev1108.py: every '; ' in the figure
# text -> ', ' (same length, so the line breaks of textwrap.fill and the layout are unchanged); output folder Figure_S2_2026-11-08. Earlier header:
# _rev1108 copy (user 2026-10-27, Figures_all_2026-11-08) of Manuscript_Figures/code/sa_tables_flowchart.py (FigS_preprocessing_flowchart.png = Figure S2, MD5).
# Flowchart part only (no CSV written). Change (REV1108_MODE=rev): 'HHCFPC (second layer)' removed from the _F_ERM output box (the only HHCFPC /
# cross-frequency text of the figure); everything else unchanged. Output FigureS2_noHHCFPC.png/.pdf (600 dpi); repro: _check/repro_FigS2.png.
"""S-A completion (user 2026-10-18, section 4.3): preprocessing flowchart, ICA statistics real vs sham, and a check of Table S2 against
the final parameters of HANDOFF section 1-1 (Revision_2026\\Preproc\\HANDOFF_20261012.md).
Inputs: Preproc\\QC\\QC_summary_wide.csv (all processed recordings), Preproc\\ica_overrides.csv (manual ICA decisions).
Outputs: data/SA_ICA_stats.csv, data/SA_parameter_check.csv, FigS_preprocessing_flowchart.png/.tiff (600 dpi)."""
data_dir = r"path/to/data"   # release: set this path
import textwrap, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import sys as _s0; _s0.path.insert(0, data_dir + r'\Revision_2026\Manuscript_Figures\code')
import tep_io as io
import os, json
from pathlib import Path
# rev1108b: own mode switch and output folder (nothing written to Figures_all_2026-11-08). REV1108B_MODE=repro: the rev1108 drawing unchanged
# (FigureS2_noHHCFPC, with semicolons) -> _check/repro_FigureS2_noHHCFPC.png; rev (default): every "; " in the figure text -> ", ".
MODE = os.environ.get("REV1108B_MODE", "rev"); assert MODE in ("rev", "repro")
S2I = os.environ.get("S2I_MODE", "rev"); assert S2I in ("rev", "start")   # 1109i
OUT = Path(data_dir + r"\Revision_2026\Manuscript_Figures\Figures_2026-11-09i"); CK = OUT / "_check"; CK.mkdir(parents=True, exist_ok=True)   # 1108n: output folder
CHANGES = []
def semi(t):   # rev1108b: semicolon -> comma + space (the rev1108 strings use "; " everywhere); records each changed line
    if MODE == "rev" and ";" in t:
        n = t.replace("; ", ", "); assert ";" not in n, t; CHANGES.append(dict(before=t, after=n)); return n
    return t

# rev1108: the ICA statistics and Table S2 check (data/SA_ICA_stats.csv, data/SA_parameter_check.csv) are not recomputed or written; flowchart only.
# ---- flowchart
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 6.5})
steps = ["Raw recording (.cnt, 64 channels incl. EOG and mastoids)", "TMS pulse detection on FC3 (threshold: real 1000, sham 300 µV/ms) and intensity labels",
         "Epochs −2 to 2 s locked to the pulse; baseline −1000 to −10 ms", "Pulse window −5 to 10 ms removed", "Bad channels (robust z > 5) and bad epochs (joint probability, 5 SD)",
         "ICA 1 (PCA 40, FastICA): decay rule + TMS-muscle classifier; manual review", "Cubic interpolation; 0.1–80 Hz band-pass, 59–61 Hz band-stop; cubic interpolation",
         "ICA 2 (FastICA): blink, eye movement, persistent muscle, electrode noise; protection rule; manual review",
         "Linear interpolation; EOG removed; channels interpolated; average reference; baseline −100 to −6 ms",
         "Condition exclusion: residual 11–25 ms > 50 µV (20–50 µV: sensitivity analysis)"]
outs = [("_F_ERM (0.1–80 Hz, 61 scalp channels)", (
         "First-layer CEEMDAN (IMFs) → ERM (IMFs 4–7), ITPC; single-participant TEPs (S-D); alpha power")),   # rev1108b: the rev1108 (noHHCFPC) text is always used
        ("_F_TEP (30 Hz low-pass)", "Conventional TEP; pre-stimulus statistics (S-B); reliable-response analysis real vs sham (S-F)")]
steps = [semi(x) for x in steps]; outs = [(semi(h), semi(d)) for h, d in outs]   # rev1108b
# rev1108c (user decision 2026-10-27, options (b)): explicit two-line splits; step 7 also "band-pass, 59–61 Hz" -> "band-pass and 59–61 Hz" (user accepted)
B_ = {"Cubic interpolation, 0.1–80 Hz band-pass, 59–61 Hz band-stop, cubic interpolation":
          "Cubic interpolation\n0.1–80 Hz band-pass and 59–61 Hz band-stop, cubic interpolation",
      "ICA 2 (FastICA): blink, eye movement, persistent muscle, electrode noise, protection rule, manual review":
          "ICA 2 (FastICA): blink, eye movement, persistent muscle, electrode noise\nprotection rule, manual review",
      "First-layer CEEMDAN (IMFs) → ERM (IMFs 4–7), ITPC, single-participant TEPs (S-D), alpha power":
          "First-layer CEEMDAN (IMFs) → ERM (IMFs 4–7), ITPC\nsingle-participant TEPs (S-D), alpha power"}
if MODE == "rev":
    for k_, v_ in B_.items(): assert sum(x == k_ for x in steps) + sum(d == k_ for _, d in outs) == 1, k_
    steps = [B_.get(x, x) for x in steps]; outs = [(h, B_.get(d, d)) for h, d in outs]
    for c_ in CHANGES: c_["after_b"] = B_.get(c_["after"], c_["after"])
# 1108n (user 2026-10-28, C): step 7 as one line without the repeated wording; output boxes renamed and their text updated (section letters S-A, S-D)
if MODE == "rev":
    _i7 = [k for k, x in enumerate(steps) if x.startswith("Cubic interpolation")]; assert len(_i7) == 1
    steps[_i7[0]] = "Cubic interpolation, 0.1–80 Hz band-pass and 59–61 Hz band-stop, cubic interpolation"
    assert [h for h, _ in outs] == ["_F_ERM (0.1–80 Hz, 61 scalp channels)", "_F_TEP (30 Hz low-pass)"], outs
    outs = [("Unfiltered version (0.1–80 Hz, 61 scalp channels)", "First-layer CEEMDAN (IMFs) → ERM (IMFs 4–7), ITPC, single-participant TEPs (S-A), alpha power"),
            ("30 Hz low-pass version", "Conventional TEP, pre-stimulus statistics (S-B), reliable-response analysis real vs sham (S-D)")]
# 1109i (user 2026-10-28): Supplementary sections numbered from 1 (S-A = S9, S-B = S10, S-D = S12)
if MODE == "rev" and S2I == "rev":
    R9 = [("single-participant TEPs (S-A)", "single-participant TEPs (S9)"), ("pre-stimulus statistics (S-B)", "pre-stimulus statistics (S10)"),
          ("reliable-response analysis real vs sham (S-D)", "reliable-response analysis real vs sham (S12)")]
    for a9, b9 in R9:
        assert sum(d.count(a9) for _, d in outs) == 1, a9
        outs = [(h, d.replace(a9, b9)) for h, d in outs]
fill_ = lambda t, w: "\n".join(textwrap.fill(x, w) for x in t.split("\n"))   # rev1108c: textwrap.fill per given line (= fill(t, w) when t has no line break)
fig, ax = plt.subplots(figsize=(6.3, 6.4)); ax.set_xlim(0, 10); ax.set_ylim(.6, len(steps) + 3.0); ax.axis("off")
y = len(steps) + 2.6
for k, s in enumerate(steps):
    yy = y - k; ax.add_patch(FancyBboxPatch((1.0, yy - .32), 8.0, .64, boxstyle="round,pad=0.02,rounding_size=0.12", fc="#F2F2F2", ec="0.35", lw=.7))
    ax.text(5.0, yy, fill_(f"{k + 1}. {s}", 92), ha="center", va="center")
    ax.annotate("", xy=(5.0, yy - .34 - .02), xytext=(5.0, yy - .66 + .02), arrowprops=dict(arrowstyle="<-", lw=.7, color="0.35"))
yb = y - len(steps) - .25
for j, (h, d) in enumerate(outs):
    x0 = .3 + j * 4.85; ax.add_patch(FancyBboxPatch((x0, yb - .75), 4.55, 1.2, boxstyle="round,pad=0.02,rounding_size=0.12", fc="white", ec="0.35", lw=.8))
    ax.text(x0 + 2.275, yb + .2, h, ha="center", va="center", weight="bold"); ax.text(x0 + 2.275, yb - .3, fill_(d, 58), ha="center", va="center", fontsize=6)
ax.text(5.0, .7, semi("Real and sham recordings: identical pipeline except the pulse-detection threshold."), ha="center", va="bottom", style="italic")
fig.tight_layout()
# rev1108b: same save procedure as rev1108_helpers.save_rev (600 dpi PNG, white, converted to RGB and tagged 600 dpi; vector PDF)
from PIL import Image
stem = ("FigureS2_1109i" if S2I == "rev" else "start_FigureS2_1108n") if MODE == "rev" else "repro_FigureS2_1108n_start"; d_ = OUT if (MODE == "rev" and S2I == "rev") else CK   # 1109i   # rev1108c: new names
fig.savefig(d_ / f"{stem}.png", dpi=600, facecolor="white"); fig.savefig(d_ / f"{stem}.pdf", facecolor="white")
im = Image.open(d_ / f"{stem}.png").convert("RGB"); im.save(d_ / f"{stem}.png", dpi=(600, 600))
fig.canvas.draw(); rr = fig.canvas.get_renderer()   # text boxes (px at 600 dpi, origin top-left) for the diff check
TB = [dict(text=t.get_text(), box=[round(v) for v in (lambda b: (b.x0 * 600 / fig.dpi, (fig.bbox.height - b.y1) * 600 / fig.dpi, b.x1 * 600 / fig.dpi, (fig.bbox.height - b.y0) * 600 / fig.dpi))(t.get_window_extent(rr))]) for t in ax.texts]
json.dump(dict(changes=CHANGES, text_boxes=TB), open(CK / f"FigureS2_text_{MODE}_{S2I}.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("saved", MODE, stem, im.size, len(CHANGES), "lines changed")
