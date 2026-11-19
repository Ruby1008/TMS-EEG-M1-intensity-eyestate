"""tiff_IJP_v3.py (user 2026-11-10, part 4, Task 1 export): Figure7_IJP_v3_rgb.png -> Figure7_IJP_v3.tif with the conversion of tiff_IJP_v2.py
(= make_final_1031.py: RGBA on white -> RGB, LANCZOS to 4252 px wide, sRGB -> RSWOP CMYK, LZW, 600 dpi);
Figure7_IJP_v3_preview.png (TIFF -> RGB, 1200 px wide, as the v2 preview); Figure7_IJP_v3_embed.png (RGB render -> 2126 px wide, LANCZOS,
optimized PNG, dpi tag 300, as Figure7_IJP_v2_embed.png)."""
data_dir = r"path/to/data"   # release: set this path
from pathlib import Path
import hashlib
from PIL import Image, ImageCms
OUT = Path(data_dir + r"\Revision_2026\Manuscript_Figures\Figures_IJP_2026-11-10")
COL = Path(r"C:\Windows\System32\spool\drivers\color")
XF = ImageCms.buildTransform(ImageCms.getOpenProfile(str(COL / "sRGB Color Space Profile.icm")), ImageCms.getOpenProfile(str(COL / "RSWOP.icm")), "RGB", "CMYK", renderingIntent=0)
XB = ImageCms.buildTransform(ImageCms.getOpenProfile(str(COL / "RSWOP.icm")), ImageCms.getOpenProfile(str(COL / "sRGB Color Space Profile.icm")), "CMYK", "RGB", renderingIntent=0)
W_PX, PPI = 4252, 600
F_TIF, F_PRE, F_EMB = OUT / "Figure7_IJP_v3.tif", OUT / "Figure7_IJP_v3_preview.png", OUT / "Figure7_IJP_v3_embed.png"
for f in (F_TIF, F_PRE, F_EMB): assert not f.exists(), f
im = Image.open(OUT / "Figure7_IJP_v3_rgb.png"); print("source", im.size, im.mode)
if im.mode in ("RGBA", "LA", "P"):
    im = im.convert("RGBA"); bg = Image.new("RGB", im.size, "white"); bg.paste(im, mask=im.split()[-1]); im = bg
else: im = im.convert("RGB")
ImageCms.applyTransform(im.resize((W_PX, round(im.size[1] * W_PX / im.size[0])), Image.LANCZOS), XF).save(F_TIF, compression="tiff_lzw", dpi=(PPI, PPI))
T = Image.open(F_TIF); T.load()
p = ImageCms.applyTransform(T, XB); p.resize((1200, round(p.size[1] * 1200 / p.size[0])), Image.LANCZOS).save(F_PRE)
im.resize((2126, round(im.size[1] * 2126 / im.size[0])), Image.LANCZOS).save(F_EMB, optimize=True, dpi=(300, 300))
print(F_TIF.name, T.size, T.mode, T.info.get("dpi"), T.info.get("compression"), f"{T.size[0] / PPI * 25.4:.1f} x {T.size[1] / PPI * 25.4:.1f} mm")
for f in (F_TIF, F_PRE, F_EMB):
    e = Image.open(f); print(f.name, e.size, e.mode, f.stat().st_size, "bytes", hashlib.sha256(f.read_bytes()).hexdigest())
