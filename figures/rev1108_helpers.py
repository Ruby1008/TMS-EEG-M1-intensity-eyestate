"""Shared helpers for the _rev1108 Supplementary figure copies (user task 2026-10-27, Figures_all_2026-11-08). Drawing helpers only; no data are read or written.
MODE: environment REV1108_MODE = "repro" (draw exactly as the original script, written to _check/repro_<stem>.png with the original savefig call, for the
pixel check against the PNG in the FINAL docx) or "rev" (default, the new marking).
  dashed_window   two black 0.8-pt dashed vertical lines at the start and end of an analysis window, label at the top of the panel above the interval
  peak_labels     component labels (text only) at the peaks, rule of the Figure 3 rev1108 copy (extreme of the curves within the component window,
                  4 pt above (P) / below (N) the extreme of the curves under the label, moved in 1-ms steps within the window if the label would touch the
                  axes frame, a red bar, the black removal bar, another label (3 pt) or come closer than 3 pt to a curve)
  signature       md5 of the x / y data of every curve, filled band and bar (for the comparison of the plotted data between the repro and rev versions)
  gray_count      number of large-area grey pixels (R = G = B, 100-245) after a 5 x 5 erosion"""
data_dir = r"path/to/data"   # release: set this path
import os, json, hashlib
from pathlib import Path
import numpy as np
import matplotlib.transforms as mtrans

MODE = os.environ.get("REV1108_MODE", "rev"); assert MODE in ("rev", "repro"), MODE
OUT = Path(data_dir + r"\Revision_2026\Manuscript_Figures\Figures_all_2026-11-08"); CK = OUT / "_check"; CK.mkdir(parents=True, exist_ok=True)

def dashed_window(ax, a, b, label, fontsize=None, ha="center", x=None, pad_pt=2, outside=False, y_data=None):
    """Black 0.8-pt dashed lines at a and b; label at the top edge of the panel, above the interval (x = centre unless given): inside the frame
    (outside=False) or just above it (outside=True, used where the curves reach the top of the panel)."""
    for v in (a, b): ax.axvline(v, color="k", lw=.8, ls=(0, (4, 2)), zorder=1, gid="rev1108_dash")
    tr = mtrans.blended_transform_factory(ax.transData, ax.transAxes) + mtrans.ScaledTranslation(0, (pad_pt if outside else -pad_pt) / 72, ax.figure.dpi_scale_trans)
    kw = {} if fontsize is None else dict(fontsize=fontsize)
    kw["bbox"] = dict(facecolor="white", edgecolor="none", pad=0.6)   # white ground so that the dashed lines do not run through the text
    if y_data is not None:   # label top pad_pt below a given height (data units), e.g. below a row of significance bars at the top of the panel
        tr = ax.transData + mtrans.ScaledTranslation(0, -pad_pt / 72, ax.figure.dpi_scale_trans)
        return ax.text((a + b) / 2 if x is None else x, y_data, label, transform=tr, ha=ha, va="top", zorder=5, gid="rev1108_winlabel", **kw)
    return ax.text((a + b) / 2 if x is None else x, 1.0, label, transform=tr, ha=ha, va="bottom" if outside else "top", zorder=5, gid="rev1108_winlabel", clip_on=False, **kw)

def _curves_px(ax, curves):
    return [ax.transData.transform(np.c_[np.interp(xx := np.linspace(t[0], t[-1], 20 * len(t)), t, t), np.interp(xx, t, y)]) for t, y in curves]

def box_curve_dist_pt(fig, ax, bb, curves):
    """Shortest distance (pt) between a display bbox and the curves (list of (t, y)), minus half of a 1.5-pt line."""
    q = _curves_px(ax, curves)
    return min(np.hypot(np.maximum(np.maximum(bb.x0 - p[:, 0], p[:, 0] - bb.x1), 0), np.maximum(np.maximum(bb.y0 - p[:, 1], p[:, 1] - bb.y1), 0)).min() for p in q) * 72 / fig.dpi - .75

def peak_labels(fig, ax, curves, windows, stars=None, fontsize=12, fontweight="bold", avoid=None, rec=None):
    """curves: list of (t, y) (same t); windows: list of (component, start_ms, end_ms); avoid: list of (x0, x1) data intervals the label must not overlap
    (e.g. the black removal bar, which spans the full height). Returns a list of dicts with the final positions."""
    stars = stars or {}; avoid = avoid or []; out = []
    fig.canvas.draw(); rr = fig.canvas.get_renderer()
    tv = curves[0][0]; G = np.array([y for _, y in curves]); lo, hi = ax.get_ylim(); x0l, x1l = ax.get_xlim()
    ptu = (hi - lo) / (ax.get_window_extent(rr).height * 72 / fig.dpi)
    red = [l.get_window_extent(rr) for l in ax.lines if l.get_linewidth() >= 3]
    placed = []
    for comp, a0, a1 in windows:
        pos = comp.startswith("P"); ww = (tv >= a0) & (tv <= a1); M = G[:, ww]; k = np.unravel_index(np.argmax(M) if pos else np.argmin(M), M.shape)
        xpk, ypk = float(tv[ww][k[1]]), float(M[k]); ci = int(k[0])
        t = ax.text(xpk, ypk, comp + stars.get(comp, ""), ha="center", va="bottom" if pos else "top", fontsize=fontsize, fontweight=fontweight, zorder=5, gid=f"lab_{comp}")
        wd = t.get_window_extent(rr).width / ax.get_window_extent(rr).width * (x1l - x0l); ok = False
        for extra in range(0, 9):
            for x_ in sorted(range(int(a0), int(a1) + 1), key=lambda v: (abs(v - xpk), v)):
                span = (tv >= x_ - wd / 2) & (tv <= x_ + wd / 2); ext = G[:, span].max() if pos else G[:, span].min(); off = (4 + extra) * (1 if pos else -1)
                t.set_position((x_, ext + off * ptu)); bb = t.get_window_extent(rr); ab = ax.get_window_extent(rr)
                bad = bb.y0 < ab.y0 + 1 or bb.y1 > ab.y1 - 1 or bb.x0 < ab.x0 + 1 or bb.x1 > ab.x1 - 1 or box_curve_dist_pt(fig, ax, bb, curves) < 3
                bad = bad or any(rb.x1 > bb.x0 and rb.x0 < bb.x1 and rb.y1 > bb.y0 - 2 and rb.y0 < bb.y1 + 2 for rb in red)
                bad = bad or any((x_ + wd / 2 > v0) and (x_ - wd / 2 < v1) for v0, v1 in avoid)
                bad = bad or any(min(bb.y1, o.y1) > max(bb.y0, o.y0) - 3 * fig.dpi / 72 and min(bb.x1, o.x1) > max(bb.x0, o.x0) - 3 * fig.dpi / 72 for o in placed)
                if not bad: ok = True; break
            if ok: break
        assert ok, comp
        placed.append(t.get_window_extent(rr))
        out.append(dict(component=comp + stars.get(comp, ""), window_ms=f"{a0}-{a1}", peak_time_ms=xpk, peak_uV=ypk, peak_curve=ci, label_time_ms=float(x_),
                        label_anchor_uV=float(ext + off * ptu), offset_pt=abs(off), gap_to_curve_pt=round(box_curve_dist_pt(fig, ax, t.get_window_extent(rr), curves), 2)))
    if rec is not None: rec.extend(out)
    return out

def label_clearance(fig, ax, text, curves):
    """Distance (pt) from a text box to the curves, and whether it overlaps another text or a full-height black bar of the panel."""
    fig.canvas.draw(); rr = fig.canvas.get_renderer(); bb = text.get_window_extent(rr)
    others = [t.get_window_extent(rr) for t in ax.texts + [ax.title] if t is not text and t.get_text()]
    ov = any(min(bb.x1, o.x1) > max(bb.x0, o.x0) and min(bb.y1, o.y1) > max(bb.y0, o.y0) for o in others)
    return dict(text=text.get_text(), gap_to_curve_pt=round(box_curve_dist_pt(fig, ax, bb, curves), 2) if curves else None, overlaps_text=bool(ov))

def _md5(a): return hashlib.md5(np.ascontiguousarray(np.asarray(a, float)).round(9).tobytes()).hexdigest()

def signature(fig):
    """Data of every plotted curve (Line2D with > 2 points), filled band (PolyCollection / FillBetween), bar (Rectangle with non-zero height) and scatter,
    per axes, without colour / line style (the marking may change, the data may not)."""
    S = []
    for ia, ax in enumerate(fig.axes):
        for l in ax.lines:
            x, y = l.get_xdata(orig=True), l.get_ydata(orig=True)
            if np.size(x) > 2 and l.get_gid() != "rev1108_dash": S.append(("line", ia, _md5(x), _md5(y)))
        for c in ax.collections:
            for p in c.get_paths():
                v = p.vertices
                if len(v) > 4: S.append(("poly", ia, _md5(v[:, 0]), _md5(v[:, 1])))
            if hasattr(c, "get_offsets") and len(c.get_offsets()) and type(c).__name__ == "PathCollection":
                o = np.asarray(c.get_offsets()); S.append(("scatter", ia, _md5(o[:, 0]), _md5(o[:, 1])))
        for p in ax.patches:
            if type(p).__name__ == "Rectangle" and p.get_height() != 0 and p.get_width() != 0 and np.isfinite(p.get_height()):
                tr = p.get_transform() - ax.transData   # data coordinates only for bars
                if tr.is_affine and np.allclose(tr.get_matrix(), np.eye(3)): S.append(("bar", ia, _md5([p.get_x(), p.get_width()]), _md5([p.get_y(), p.get_height()])))
    return sorted(S)

def save_signature(fig, stem):
    json.dump(signature(fig), open(CK / f"sig_{MODE}_{stem}.json", "w"))

def gray_count(png):
    from PIL import Image
    from scipy.ndimage import binary_erosion
    a = np.asarray(Image.open(png).convert("RGB")).astype(int)
    g = (a[..., 0] == a[..., 1]) & (a[..., 1] == a[..., 2]) & (a[..., 0] >= 100) & (a[..., 0] <= 245)
    return int(binary_erosion(g, structure=np.ones((5, 5), bool)).sum()), int(g.sum())

def save_rev(fig, stem, **kw):
    """600 dpi PNG (RGB, white) + vector PDF in the output folder."""
    from PIL import Image
    fig.savefig(OUT / f"{stem}.png", dpi=600, facecolor="white", **kw); fig.savefig(OUT / f"{stem}.pdf", facecolor="white", **kw)
    im = Image.open(OUT / f"{stem}.png").convert("RGB"); im.save(OUT / f"{stem}.png", dpi=(600, 600))

def gray_without(fig, pred, stem, **kw):
    """Diagnostic render to _check/<stem>_diag.png with the lines for which pred(line) is True hidden (e.g. the grey single-channel lines), so that the grey
    count can be attributed; the lines are made visible again afterwards."""
    hid = [l for ax in fig.axes for l in ax.lines if pred(l)]
    for l in hid: l.set_visible(False)
    p = CK / f"{stem}_diag.png"; fig.savefig(p, dpi=600, facecolor="white", **kw)
    for l in hid: l.set_visible(True)
    return len(hid), gray_count(p)
