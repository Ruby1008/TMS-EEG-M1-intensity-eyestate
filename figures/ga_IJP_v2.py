"""ga_IJP_v2.py (user 2026-11-10, part 6, Task A): copy of ga_IJP.py with TEXT changes only. Plotted data, colours, marker/line sizes and canvas
unchanged. Changed: title, panel titles 2 and 3, panel-1 label "motor threshold" -> "RMT", RMT label and "no sham MEP recorded" 24 -> 28 px,
two new bottom lines + one grey 22-px note line; to fit the third bottom line the bottom band grows 120 -> 160 px (panel height -40 px, fonts unchanged).
Outputs go to Desktop part6 folder with the _v2 suffix. Original docstring follows.
ga_IJP.py (user 2026-11-10, part 5, Task B): graphical abstract, three readouts of the M1 TMS pulse across intensity, real vs sham.
No new statistics: group mean +/- SE (SD / sqrt(n), ddof = 1) of existing per-participant values, eye states averaged within participant.
  Panel 1 MEP: MEP_reprocess_2026-10-22/main_qc2/MEP_matrix_reprocessed.csv, columns EO_80..EC_120; ln of each cell, then mean of EO and EC
          (as logmep_tables_1112.py); real only, n = 23.
  Panel 2 early EEG response: data/R21_r15_cells.csv, column log_snr_local, 17 main participants (io.MAIN17) = the values drawn in Figure 4B
          (fig4_IJP.py). log SNR = ln(GMFP 15-60 ms / mean GMFP of the null), so 0 = no response. Participants with both eye states only
          (real 17; sham 16, 024 has sham EC only). GA_SNR_SOURCE=r21 switches to data/R21_cells.csv, column log_snr (original null).
  Panel 3 alpha phase locking: Review_2026-11-10/N2/N2_sham_ITPC_summary.csv, rows IMF = 6 (window 27-456 ms), column ITPC_mean; real = all 17,
          sham = rows with in_16 = True (16; 024 has sham EC only).
Canvas 2656 x 1062 px (13.28 x 5.31 in at 200 dpi), Arial. Outputs (Desktop part5 folder): GA_IJP.png (RGB), GA_IJP.tif (same pixels, 300 dpi tag,
RGB, LZW), GA_IJP.pdf (vector), GA_IJP_notext.png (bands empty), GA_IJP_values.csv, GA_IJP_check_1328.png (half size, legibility check)."""
data_dir = r"path/to/data"   # release: set this path
import os, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
MF = Path(data_dir + r"\Revision_2026\Manuscript_Figures"); sys.path.insert(0, str(MF / "code")); import tep_io as io
OUT = Path(data_dir + r"\outputs\IJP_upload_20261110_part6_drag_to_Drive")
DPI = 200; W, H = 2656, 1062; PX = 72 / DPI                      # points per pixel
REAL, SHAM = "#C5550A", "0.45"
LEV = [80, 100, 120]

def summ(df, val, panel, cond, src, col):
    rows = []
    for lv in LEV:
        x = df[df.intensity == lv][val].dropna(); rows.append(dict(panel=panel, condition=cond, intensity=lv, mean=x.mean(), SE=x.std(ddof=1) / np.sqrt(len(x)), n=len(x), source=src, column=col))
    return rows

V = []
# ---- panel 1: MEP (real)
M = pd.read_csv(MF / "MEP_reprocess_2026-10-22" / "main_qc2" / "MEP_matrix_reprocessed.csv"); assert len(M) == 23 and not M.isna().any().any()
L = pd.DataFrame({"sub": M.subject})
for lv in LEV: L[lv] = np.log(M[[f"EO_{lv}", f"EC_{lv}"]]).mean(1)
P1 = L.melt(id_vars="sub", var_name="intensity", value_name="logMEP"); P1["intensity"] = P1.intensity.astype(int)
V += summ(P1, "logMEP", "1 MEP", "Real", "MEP_reprocess_2026-10-22/main_qc2/MEP_matrix_reprocessed.csv", "EO_/EC_80..120, ln then mean of EO and EC")
# ---- panel 2: early EEG response
src2, col2 = (("data/R21_cells.csv", "log_snr") if os.environ.get("GA_SNR_SOURCE", "fig4b") == "r21" else ("data/R21_r15_cells.csv", "log_snr_local"))
C = pd.read_csv(MF / src2); C = C[C["sub"].isin(io.MAIN17)]
for st in ("Real", "Sham"):
    c = C[C.stim == st]; both = c.groupby(["sub", "intensity"]).eye.transform("nunique") == 2   # participants x intensity with EO and EC
    p = c[both].groupby(["sub", "intensity"])[col2].mean().reset_index()
    V += summ(p, col2, "2 early EEG", st, src2, col2 + ", eye states averaged within participant (both eye states required)")
# ---- panel 3: ITPC IMF6
N = pd.read_csv(MF / "Review_2026-11-10" / "N2" / "N2_sham_ITPC_summary.csv"); N = N[N.IMF == 6]; assert set(N.window_ms) == {"27-456"}
for st, sel in (("Real", N.stim == "real"), ("Sham", (N.stim == "sham") & N.in_16)):
    p = N[sel].groupby(["participant", "intensity"]).ITPC_mean.mean().reset_index()
    V += summ(p, "ITPC_mean", "3 ITPC", st, "Review_2026-11-10/N2/N2_sham_ITPC_summary.csv", "ITPC_mean, IMF 6, 27-456 ms" + (", in_16 = True" if st == "Sham" else ""))
VD = pd.DataFrame(V)
assert list(VD[VD.panel == "1 MEP"].n) == [23] * 3 and list(VD[(VD.panel == "3 ITPC") & (VD.condition == "Real")].n) == [17] * 3 and list(VD[(VD.panel == "3 ITPC") & (VD.condition == "Sham")].n) == [16] * 3

# ---- drawing
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "axes.linewidth": 1.4 * 200 / 144, "xtick.major.width": 1.2, "ytick.major.width": 1.2})
TICK, PTITLE, TITLE, BOT, NOTE, SMALL = 28 * PX, 34 * PX, 46 * PX, 36 * PX, 28 * PX, 22 * PX   # v2: NOTE 24 -> 28 px, SMALL = new grey line
def draw(text=True):
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI, facecolor="white")
    top_band, bot_band = 150, 160                                  # px
    ax_bottom_px, ax_top_px = bot_band + 150, H - top_band - 62     # room for ticks, x label, note / panel title
    lefts_px = [150, 1040, 1930]; ax_w_px = 640
    axes = [fig.add_axes([l / W, ax_bottom_px / H, ax_w_px / W, (ax_top_px - ax_bottom_px) / H]) for l in lefts_px]
    spec = (("1 MEP", "Corticospinal output (MEP)", "log MEP (ln mV)"), ("2 early EEG", "Early EEG response (15–60 ms, Fig. 4)", "log SNR"),
            ("3 ITPC", "Alpha phase locking (IMF6 ITPC)", "ITPC"))
    for k, (ax, (pan, ttl, yl)) in enumerate(zip(axes, spec)):
        ax.axvline(100, color="0.6", lw=1.6, ls=(0, (4, 3)), zorder=0)
        if pan == "2 early EEG": ax.axhline(0, color="0.6", lw=1.2, zorder=0)
        for cond, col, ls, mfc, dx in (("Real", REAL, "-", REAL, -1.6), ("Sham", SHAM, "--", "white", 1.6)):
            d = VD[(VD.panel == pan) & (VD.condition == cond)].sort_values("intensity")
            if len(d): ax.errorbar(d.intensity + dx, d["mean"], yerr=d.SE, color=col, ls=ls, lw=2.6, marker="o", ms=11, mfc=mfc, mec=col, mew=2.2, capsize=5, capthick=2, elinewidth=2,
                                   label=f"{cond} TMS")
        ax.set_xlim(70, 130); ax.set_xticks(LEV); ax.set_xticklabels(["80%", "100%", "120%"]); ax.tick_params(labelsize=TICK, length=6)
        ax.set_xlabel("Intensity (% RMT)", fontsize=TICK, labelpad=4); ax.set_ylabel(yl, fontsize=TICK, labelpad=6)
        ax.set_title(ttl, fontsize=PTITLE, pad=14); ax.spines[["top", "right"]].set_visible(False)
        if k == 0:
            ax.text(101.2, 0.98, "RMT", transform=ax.get_xaxis_transform(), fontsize=NOTE, color="0.4", va="top", ha="left")
            ax.text((lefts_px[0] + ax_w_px / 2) / W, (bot_band + 14) / H, "no sham MEP recorded", transform=fig.transFigure, fontsize=NOTE, color="0.45", ha="center", va="bottom")
    from matplotlib.lines import Line2D
    hh = [Line2D([], [], color=REAL, ls="-", lw=2.6, marker="o", ms=11, mfc=REAL, mec=REAL, mew=2.2), Line2D([], [], color=SHAM, ls="--", lw=2.6, marker="o", ms=11, mfc="white", mec=SHAM, mew=2.2)]
    axes[0].legend(hh, ["Real TMS", "Sham TMS"], loc="upper left", frameon=False, fontsize=TICK, handlelength=2.6, borderaxespad=0.2)
    if text:
        fig.text(0.5, 1 - 75 / H, "Real vs sham M1 TMS at 80–120% RMT: three readouts of one pulse", ha="center", va="center", fontsize=TITLE, fontweight="bold")
        fig.text(0.5, 120 / H, "Real, but rarely sham, TMS evoked an early EEG response at 80% RMT", ha="center", va="center", fontsize=BOT, fontweight="bold", color="0.2")
        fig.text(0.5, 74 / H, "Single-trial P30 covaried with the MEP at 100–120% RMT. No consistent effect of eye closure or mu power.", ha="center", va="center", fontsize=BOT, color="0.2")
        fig.text(0.5, 28 / H, "Mean ± SE, eye states averaged. MEP n = 23. EEG n = 17 (real), 16 (sham).", ha="center", va="center", fontsize=SMALL, color="0.45")
    return fig, axes

fig, axes = draw(True)
# overlap check on the drawn Agg canvas (before the PDF save): every drawn text element vs every other, in pixels; tick labels only for ticks inside the view
fig.canvas.draw(); r = fig.canvas.get_renderer()
tl = set(); [tl.update(id(t) for t in (a.get_xticklabels() + a.get_yticklabels())) for a in axes]
def _drawn(t, a=None):
    if id(t) not in tl: return True
    for a in axes:
        for lab, lim, loc in ((a.get_xticklabels(), a.get_xlim(), a.get_xticks()), (a.get_yticklabels(), a.get_ylim(), a.get_yticks())):
            for tt, v in zip(lab, loc):
                if tt is t: return min(lim) <= v <= max(lim)
    return True
texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_text() and t.get_visible() and _drawn(t)]
bbs = [(t.get_text(), t.get_window_extent(r)) for t in texts]
over = [(a, b) for i, (a, ba) in enumerate(bbs) for (b, bb) in bbs[i + 1:] if ba.overlaps(bb) and not {a, b} == {"Real TMS", "Sham TMS"}]   # legend rows touch by design
inside = all(0 <= bb.x0 and bb.x1 <= W and 0 <= bb.y0 and bb.y1 <= H for _, bb in bbs)
fig.savefig(OUT / "GA_IJP_v2_rgba.png", dpi=DPI, facecolor="white"); fig.savefig(OUT / "GA_IJP_v2.pdf", facecolor="white"); plt.close(fig)
im = Image.open(OUT / "GA_IJP_v2_rgba.png"); assert im.size == (W, H), im.size
im.convert("RGB").save(OUT / "GA_IJP_v2.png"); im.convert("RGB").save(OUT / "GA_IJP_v2.tif", compression="tiff_lzw", dpi=(300, 300)); os.remove(OUT / "GA_IJP_v2_rgba.png")
Image.open(OUT / "GA_IJP_v2.png").resize((1328, 531), Image.LANCZOS).save(OUT / "GA_IJP_v2_check_1328.png")
fig2, _ = draw(False); fig2.savefig(OUT / "GA_IJP_v2_notext_rgba.png", dpi=DPI, facecolor="white"); plt.close(fig2)
Image.open(OUT / "GA_IJP_v2_notext_rgba.png").convert("RGB").save(OUT / "GA_IJP_v2_notext.png"); os.remove(OUT / "GA_IJP_v2_notext_rgba.png")
VD.to_csv(OUT / "GA_IJP_v2_values.csv", index=False)
print(VD[["panel", "condition", "intensity", "mean", "SE", "n"]].round(4).to_string(index=False))
print("text overlaps:", over if over else "none", "| all text inside canvas:", inside)
for f in ("GA_IJP_v2.png", "GA_IJP_v2.tif", "GA_IJP_v2_notext.png"):
    x = Image.open(OUT / f); print(f, x.size, x.mode, x.info.get("dpi"), x.info.get("compression"))
