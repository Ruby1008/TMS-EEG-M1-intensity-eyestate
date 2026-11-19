# _1114 copy (user 2026-11-14, Figures_2026-11-14_minor) of Figures_2026-11-08v/code/fig1_1108v.py. M1114=repro: unchanged drawing -> _check/repro; M1114=rev:
# 'MEP Amplitude (mV)' -> 'MEP amplitude (mV)' (y label; font size still measured on the original string) and 'Amplitude' -> 'amplitude' (bracket label). Lines marked # 1114.
"""Main Figure 1 (experimental procedure), redrawn by code (user 2026-10-28, Figures_2026-11-08v).
Layout, positions and colours follow the figure in spTMS-EEG_Manuscript_CerCor_FINAL_2026-11-07d.docx (word/media/image1.png, 1753 x 1017 px, saved as
_check/Figure1_reference.png). The whole figure is one axes in reference pixel units (x 0..1753, y 0..1017 downwards); the canvas is 17.53 x 10.17 in at
300 dpi (5259 x 3051 px, 3x the reference), so 1 reference px = 0.72 pt.
- Text: all text is typeset here (Times New Roman); the new wording is set in TEXT below. Font sizes and baselines are measured from the reference lines
  (ink width of the old string in Times New Roman, _check/Figure1_text_measure.json, from measure_fig1_text_1108v.py; baseline = lowest row whose ink
  count is >= 25% of the line maximum).
- Head, EEG cap, figure-8 coil, hand and lightning bolts: vector polygons traced from the reference colour masks (filled contours of the masks, upsampled
  4x with light smoothing); the label text next to the hand is masked out before tracing. Red dots are circles (centroid and area from the reference).
  Screens, crosses, timelines, ticks and arrows are drawn as vector shapes at the measured positions (timelines drawn horizontal).
- MEP inset: the original panel waveform is trial 53 of participant 002, eyes-open block, 100% RMT (row_id 431 of
  MEP_reprocess_2026-10-22/extract/mep_waveforms.csv, file TMS_MEP_002_1230_EO_acq_to_Veta13-Apr-2023_preprocessed.mat; r = 0.99999999 with the
  waveform recovered from Figure1_MEP_panel_editable.pptx, find_fig1_mep_trial_1108v.py). The extract is trigger-aligned and the TMS pulse is at -11 ms
  (pulse_ms_re_trigger in main_qc2/mep_trials_reprocessed.csv), so time is shifted by +11 ms (0 ms = TMS pulse), as in the original panel.
Outputs: Figure1_1108v.png, Figure1_1108v.pdf."""
data_dir = r"path/to/data"   # release: set this path
import json
from pathlib import Path
import numpy as np, pandas as pd
from PIL import Image
from scipy import ndimage
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Polygon

import os as _os; M1114 = _os.environ.get("M1114", "rev"); assert M1114 in ("repro", "rev")   # 1114
OUT = Path(data_dir + r"\Revision_2026\Manuscript_Figures\Figures_2026-11-08v"); CK = OUT / "_check"   # inputs (read only)
OUT14 = Path(data_dir + r"\Revision_2026\Manuscript_Figures\Figures_2026-11-14_minor"); CK14 = OUT14 / "_check" / ("repro" if M1114 == "repro" else ""); CK14.mkdir(parents=True, exist_ok=True)   # 1114: outputs
EX = Path(data_dir + r"\Revision_2026\Manuscript_Figures\MEP_reprocess_2026-10-22\extract")
plt.rcParams.update({"font.family": "Times New Roman", "pdf.fonttype": 42, "mathtext.fontset": "custom", "axes.unicode_minus": True})
PT = 0.72                                  # pt per reference px
RB, TMSC, MAG, GOLD, REDD, BLUE, COIL = "#942911", "#9b1c1c", "#c91d67", "#ffda0a", "#ff0000", "#1f77b4", "#474747"   # measured reference colours

A = np.asarray(Image.open(CK / "Figure1_reference.png").convert("RGBA")).astype(float)
REF = np.where(A[..., 3:4] > 0, A[..., :3] * A[..., 3:4] / 255 + 255 * (1 - A[..., 3:4] / 255), 255)   # composite on white
R, G, B = REF[..., 0], REF[..., 1], REF[..., 2]; MX = REF.max(2); MN = REF.min(2)
MEAS = json.load(open(CK / "Figure1_text_measure.json"))

fig = plt.figure(figsize=(17.53, 10.17), dpi=300)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1753); ax.set_ylim(1017, 0); ax.axis("off")


# ---------------------------------------------------------------- traced illustration
def trace(field, box, color, z, up=4, sigma=0.5, min_px=12):
    """Filled contour (level 0.5) of a coverage field inside box (x0, y0, x1, y1), drawn as vector polygons."""
    x0, y0, x1, y1 = box
    f = field[y0:y1, x0:x1].astype(float)
    lab, n = ndimage.label(f > 0.5)
    small = np.isin(lab, [i for i, c in enumerate(ndimage.sum(np.ones_like(f), lab, range(n + 1))) if i > 0 and c < min_px])
    f[small] = 0
    f = np.pad(ndimage.gaussian_filter(f, sigma), 2)
    f = ndimage.zoom(f, up, order=1)
    xs = np.linspace(x0 - 2 + 0.5, x1 + 2 - 0.5, f.shape[1]); ys = np.linspace(y0 - 2 + 0.5, y1 + 2 - 0.5, f.shape[0])
    ax.contourf(xs, ys, f, levels=[0.5, 10], colors=[color], antialiased=True, zorder=z)


ILL = (1180, 40, 1660, 380)                 # head, coil and hand
COILBOX = (1242, 40, 1346, 140)
TXTBOX = (1440, 85, 1665, 158)              # old hand label (masked out, retyped below)
black = np.clip(1 - MX / 255, 0, 1)
cx0, cy0, cx1, cy1 = COILBOX
black[cy0:cy1, cx0:cx1] = np.clip(1 - MX[cy0:cy1, cx0:cx1] / 110, 0, 1)        # inside the coil box the grey coil must not count as black
tx0, ty0, tx1, ty1 = TXTBOX; black[ty0:ty1, tx0:tx1] = 0
magenta = ((R > 150) & (G < 110) & (B > 60) & (R - G > 80)).astype(float)
mag_d = ndimage.binary_dilation(magenta > 0, iterations=1)
in_coil = np.zeros_like(magenta, bool); in_coil[cy0:cy1, cx0:cx1] = True
magenta = np.where(in_coil, mag_d, magenta > 0).astype(float)                    # extend the cap 1 px under the coil (no seam)
coil = ((MX - MN < 25) & (MX > 30) & (MX < 175)); coil[~in_coil] = False
coil = ndimage.binary_opening(coil, np.ones((3, 3))).astype(float)              # drop the anti-aliased fringe of the black outline
gold = ndimage.binary_dilation((R > 200) & (G > 170) & (B < 120), iterations=1).astype(float)

trace(magenta, (1285, 75, 1420, 165), MAG, 2)
trace(black, ILL, "black", 3)
trace(coil, COILBOX, COIL, 4, sigma=1.0)                                     # binary mask: smooth more
BOLTS = (900, 690, 1350, 812)
trace(gold, BOLTS, GOLD, 2)
trace(black, BOLTS, "black", 3)
red = (R > 200) & (G < 60) & (B < 60)
lab, n = ndimage.label(red)
for i in range(1, n + 1):
    yy, xx = np.where(lab == i)
    if len(xx) > 30 and xx.mean() > 1180:
        ax.add_patch(Circle((xx.mean() + 0.5, yy.mean() + 0.5), np.sqrt(len(xx) / np.pi) + 0.3, color=REDD, lw=0, zorder=6))

# ---------------------------------------------------------------- screens, timelines, arrows
for x0, y0, x1, y1 in [(146, 105, 530, 374), (714, 106, 1098, 374)]:
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, color="black", lw=0, zorder=1))
    cx, cy = (x0 + x1) / 2 - 0.5, 241.5
    ax.plot([cx - 17.5, cx + 17.5], [cy, cy], color="white", lw=4.5 * PT, solid_capstyle="butt", zorder=2)
    ax.plot([cx, cx], [cy - 17, cy + 17], color="white", lw=4.5 * PT, solid_capstyle="butt", zorder=2)


def arrow_line(x0, x1, y, lw_px, color, head_l=20, head_w=9, both=False, z=3):
    xa = x0 + (head_l - 2 if both else 0); xb = x1 - head_l + 2
    ax.plot([xa, xb], [y, y], color=color, lw=lw_px * PT, solid_capstyle="butt", zorder=z)
    ax.add_patch(Polygon([(x1, y), (x1 - head_l, y - head_w), (x1 - head_l, y + head_w)], closed=True, color=color, lw=0, zorder=z))
    if both:
        ax.add_patch(Polygon([(x0, y), (x0 + head_l, y - head_w), (x0 + head_l, y + head_w)], closed=True, color=color, lw=0, zorder=z))


arrow_line(68, 1583, 400, 6, "black", head_l=20, head_w=8)                     # session timeline
ax.plot([643, 643], [367, 441], color="black", lw=6 * PT, solid_capstyle="butt", zorder=3)
arrow_line(821, 1493, 847, 5.5, "black", head_l=18, head_w=7.5)                # trial timeline
for x in (969.5, 1148.5, 1325.5):
    ax.plot([x, x], [816, 879], color="black", lw=6 * PT, solid_capstyle="butt", zorder=4)
arrow_line(973, 1146, 863, 4, RB, head_l=17, head_w=8, both=True)              # ITI arrows
arrow_line(1152, 1323, 863, 4, RB, head_l=17, head_w=8, both=True)


# ---------------------------------------------------------------- text
def baseline(name):
    x0, y0, x1, y1 = MEAS[name]["ink_box"]; m = (MX[y0:y1 + 1, x0:x1 + 1] < 140) | ((R - G)[y0:y1 + 1, x0:x1 + 1] > 60)
    cnt = m.sum(1); return y0 + np.where(cnt >= 0.25 * cnt.max())[0].max() + 1


def put(name, s, x, ha, color="black", weight="normal", size=None, dy=0):
    fs = size or MEAS[name]["fontsize_pt"]
    ax.text(x, baseline(name) + dy, s, ha=ha, va="baseline", fontsize=fs, fontweight=weight, color=color, zorder=8)


def cx_(name): b = MEAS[name]["ink_box"]; return (b[0] + b[2] + 1) / 2


def lx_(name): return MEAS[name]["ink_box"][0]


TEXT = {   # name: (new string, alignment, colour, weight)
    "title1": ("Localizing M1 Hotspot and RMT", "center", "black", "bold"),
    "title2": ("Single-Pulse TMS-EEG", "center", "black", "bold"),
    "fdi1": ("First dorsal", "center", "black", "normal"),
    "fdi2": ("interosseous (FDI)", "center", "black", "normal"),
    "rmt_title": ("Resting Motor Threshold (RMT) Criteria", "center", RB, "bold"),
    "rmt1": ("At least 5 of 10 MEPs > 0.05 mV,", "center", "black", "normal"),
    "rmt2": ("stimuli delivered at least 5 s apart", "center", "black", "normal"),
    "eye1": ("Eyes closed (EC) and eyes open (EO) blocks", "left", "black", "normal"),
    "eye2": ("15 min each, counterbalanced order", "left", "black", "normal"),
    "int_title": ("TMS intensity:", "left", "black", "bold"),
    "int_line": ("80%, 100% and 120% RMT (interleaved)", "left", "black", "normal"),
    "trials": ("126 trials per block", "left", "black", "normal"),
    "iti": ("ITI", "center", RB, "bold"),
    "iti_val": ("5\u20137 s", "center", RB, "normal"),
    "tms": ("TMS", "center", TMSC, "bold"),
    "mep1": ("MEP", "center", "black", "normal"),
    "mep2": ("Amplitude" if M1114 == "repro" else "amplitude", "center", "black", "normal"),   # 1114
    "xlabel": ("Time (ms)", "center", "black", "bold"),
}
for k, (s, ha, col, wt) in TEXT.items():
    x = cx_(k) if ha == "center" else lx_(k)
    if k == "tms": x = 179.5                                                     # centred on the 0-ms line
    put(k, s, x, ha, col, wt)
# y label (rotated): same size as in the reference (ink height 636..887 of "MEP Amplitude (mV)", bold)
f2 = plt.figure(figsize=(8, 1), dpi=400); f2.text(0.01, 0.3, "MEP Amplitude (mV)", fontsize=10, fontweight="bold"); f2.canvas.draw()
im = np.asarray(f2.canvas.buffer_rgba())[..., 0]; cc = np.where((im < 128).any(0))[0]; plt.close(f2)
ys_fs = 10 * ((887 - 636 + 1) / 100) / ((cc.max() - cc.min() + 1) / 400)
ax.text(28.5, (636 + 888) / 2, "MEP Amplitude (mV)" if M1114 == "repro" else "MEP amplitude (mV)", rotation=90,   # 1114
        ha="center", va="center", fontsize=ys_fs, fontweight="bold", zorder=8)

# ---------------------------------------------------------------- MEP inset
W = pd.read_csv(EX / "mep_waveforms.csv"); row = W[W.row_id == 431].iloc[0]
SHIFT = 11
t_ex = np.arange(-12 - SHIFT, 72 - SHIFT + 1); y = row[[f"t{t:+04d}" for t in t_ex]].to_numpy(float); t = t_ex + SHIFT
X0, X1, PXMS = -12, 72, (618.5 - 179.5) / 70; Y0, Y1, PXMV = -0.72, 0.72, (915 - 625) / 1.2
left = 179.5 + X0 * PXMS; right = 179.5 + X1 * PXMS; top = 770 - Y1 * PXMV; bot = 770 - Y0 * PXMV
ia = fig.add_axes([left / 1753, 1 - bot / 1017, (right - left) / 1753, (bot - top) / 1017])
ia.set_xlim(X0, X1); ia.set_ylim(Y0, Y1); ia.patch.set_alpha(0)
for s in ("top", "right"): ia.spines[s].set_visible(False)
for s in ("left", "bottom"): ia.spines[s].set_linewidth(3 * PT)
ia.plot(t, y, color=BLUE, lw=3.2 * PT, zorder=3, solid_joinstyle="round")
ia.axvline(0, color=TMSC, lw=4 * PT, dashes=(13 / 4, 12 / 4), zorder=2)
ymax, ymin = y.max(), y.min(); xb = (387 - 179.5) / PXMS
ia.plot([xb - 0.9, xb, xb, xb - 0.9], [ymax, ymax, ymin, ymin], color="black", lw=4 * PT, solid_joinstyle="miter", solid_capstyle="butt", zorder=4)
xt = np.arange(-10, 71, 10); yt = np.round(np.arange(-0.6, 0.61, 0.2), 1)
ia.set_xticks(xt); ia.set_xticklabels([f"{v:d}".replace("-", "\u2212") for v in xt])
ia.set_yticks(yt); ia.set_yticklabels([f"{v:.1f}".replace("-", "\u2212") if v != 0 else "0.0" for v in yt])
tfs = (MEAS["ytick"]["fontsize_pt"] + MEAS["xtick"]["fontsize_pt"]) / 2
ia.tick_params(direction="out", length=3.5 * PT, width=2 * PT, labelsize=tfs)
ia.tick_params(axis="y", pad=(104 - 4 - 94) * PT); ia.tick_params(axis="x", pad=(952 - 949) * PT)

_stem = CK14 / "repro_Figure1_1108v" if M1114 == "repro" else OUT14 / "Figure1_1114"   # 1114
fig.savefig(f"{_stem}.png", dpi=300)
fig.savefig(f"{_stem}.pdf")
json.dump(dict(waveform_row_id=431, shift_ms=SHIFT, t_ms=[int(v) for v in t], mV=[float(v) for v in y], p2p=float(ymax - ymin), ylabel_fs=float(ys_fs),
               tick_fs=float(tfs)), open(OUT14 / "_check" / f"Figure1_mep_inset_{M1114}.json", "w"), indent=1)
print("saved; p2p", round(ymax - ymin, 4), "ylabel fs", round(ys_fs, 2), "tick fs", round(tfs, 2))
