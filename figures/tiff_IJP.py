"""tiff_IJP.py (user 2026-11-10, Part F3/F4): Figure{4,7}_IJP_rgb.png -> CMYK TIFF with the conversion of make_final_1031.py
(RGBA on white -> RGB, LANCZOS resize to 4252 px wide = 180 mm at 600 ppi, sRGB -> RSWOP ICC transform, LZW, dpi 600) -> Figure{4,7}_IJP.tif;
previews Figure{4,7}_IJP_preview.png (TIFF -> RGB with the inverse transform, LANCZOS to 1600 px wide). Also reports the TIFF specification
and, as a check of the conversion, converts the current PNGs (Figure4_1111, Figure7_1114) the same way in memory: outside the changed title box
the IJP TIFF must equal that conversion (0 differing pixels)."""
data_dir = r"path/to/data"   # release: set this path
from pathlib import Path
import numpy as np, pandas as pd
from PIL import Image, ImageCms

MF = Path(data_dir + r"\Revision_2026\Manuscript_Figures"); OUT = MF / "Figures_IJP_2026-11-10"
COL = Path(r"C:\Windows\System32\spool\drivers\color")
XF = ImageCms.buildTransform(ImageCms.getOpenProfile(str(COL / "sRGB Color Space Profile.icm")), ImageCms.getOpenProfile(str(COL / "RSWOP.icm")), "RGB", "CMYK", renderingIntent=0)
XB = ImageCms.buildTransform(ImageCms.getOpenProfile(str(COL / "RSWOP.icm")), ImageCms.getOpenProfile(str(COL / "sRGB Color Space Profile.icm")), "CMYK", "RGB", renderingIntent=0)
W_MM, PPI = 180, 600; W_PX = round(W_MM / 25.4 * PPI)
def to_cmyk(p):   # = make_final_1031.py
    im = Image.open(p)
    if im.mode in ("RGBA", "LA", "P"):
        im = im.convert("RGBA"); bg = Image.new("RGB", im.size, "white"); bg.paste(im, mask=im.split()[-1]); im = bg
    else: im = im.convert("RGB")
    h_px = round(im.size[1] * W_PX / im.size[0]); im = im.resize((W_PX, h_px), Image.LANCZOS)
    return ImageCms.applyTransform(im, XF)
rows = []
for name, ref in (("Figure4", MF / "Figure4_2026-11-11" / "Figure4_1111.png"), ("Figure7", MF / "Figures_2026-11-14_minor" / "Figure7_1114.png")):
    c = to_cmyk(OUT / f"{name}_IJP_rgb.png"); f = OUT / f"{name}_IJP.tif"; assert not f.exists(), f
    c.save(f, compression="tiff_lzw", dpi=(PPI, PPI))
    T = Image.open(f); T.load()
    p = ImageCms.applyTransform(T, XB); p.resize((1600, round(p.size[1] * 1600 / p.size[0])), Image.LANCZOS).save(OUT / f"{name}_IJP_preview.png")
    a = np.asarray(T).astype(int); b = np.asarray(to_cmyk(ref)).astype(int); m = np.abs(a - b).sum(2) > 0; ys, xs = np.where(m)
    rows.append(dict(file=f.name, width_px=T.size[0], height_px=T.size[1], mode=T.mode, dpi=str(T.info.get("dpi")), compression=T.info.get("compression"),
                     width_mm=round(T.size[0] / PPI * 25.4, 1), height_mm=round(T.size[1] / PPI * 25.4, 1), file_bytes=f.stat().st_size,
                     diff_px_vs_current_png=int(m.sum()), diff_box=f"{xs.min()},{ys.min()},{xs.max()},{ys.max()}" if m.any() else "none",
                     preview=f"{name}_IJP_preview.png", preview_px="x".join(map(str, Image.open(OUT / f'{name}_IJP_preview.png').size))))
R = pd.DataFrame(rows); R.to_csv(OUT / "_check" / "tiff_spec_IJP.csv", index=False); print(R.to_string())
# Figure 7 panel C: horizontal gaps between the new title line and the panel letter 'C' / the colorbar label 't' (dark pixels, RGB render, 600 dpi)
G = np.asarray(Image.open(OUT / "Figure7_IJP_rgb.png").convert("L")) < 128
band = G[3086:3155]; cols = np.where(band.any(0))[0]; title = cols[(cols >= 150) & (cols <= 1250)]
letter = np.where(G[3100:3300, 0:160].any(0))[0]; tlab = np.where(G[3150:3300, 1250:1400].any(0))[0] + 1250
print("Figure7 title line x", title.min(), "-", title.max(), "| letter C right edge x", letter.max(), "-> gap", title.min() - letter.max() - 1, "px",
      "| colorbar 't' left edge x", tlab.min() if len(tlab) else None)
