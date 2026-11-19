"""Shared helpers (2026-10-19): the 2 x 3 rmANOVA summary, BH-FDR and run formatting used by tep_stats.py (copied so that the new scripts do not re-run it)."""
import numpy as np, pandas as pd, pingouin as pg

def bh(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); q = np.empty(m)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]; return np.minimum(q, 1)

def anova2(df, dv):
    a = pg.rm_anova(data=df, dv=dv, within=["eye", "intensity"], subject="sub", correction=True, effsize="np2").set_index("Source"); out = {}
    for src, key in (("eye", "eye"), ("intensity", "intensity"), ("eye * intensity", "interaction")):
        r = a.loc[src]; pgg = r["p_GG_corr"] if "p_GG_corr" in r and pd.notna(r.get("p_GG_corr")) else r["p_unc"]
        out.update({f"F_{key}": r["F"], f"df1_{key}": r["ddof1"], f"df2_{key}": r["ddof2"], f"p_{key}": pgg, f"np2_{key}": r["np2"], f"eps_{key}": r["eps"] if "eps" in r and pd.notna(r.get("eps")) else 1.0})
    return out

def runs(starts, mask, width=20):
    out, st, prev = [], None, None
    for a0, m in zip(starts, mask):
        if m and st is None: st = a0
        if not m and st is not None: out.append(f"{st}-{prev + width}"); st = None
        prev = a0
    if st is not None: out.append(f"{st}-{prev + width}")
    return ", ".join(out) or "none"
