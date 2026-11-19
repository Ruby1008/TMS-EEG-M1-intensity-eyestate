"""Shared helpers for Review_2026-11-10 (user 2026-10-28). Reads the main-analysis inputs exactly as Manuscript_Figures/code/tep_stats.py does
(17 participants, ROI FC3/FC1/C3/C1, ERM = erm_<id>.mat erm_roi, conventional TEP = data/TEP_roi_waveforms.csv file 'TEP' = _F_TEP 0.1-30 Hz),
plus APA formatting. No input file is modified."""
data_dir = r"path/to/data"   # release: set this path
import sys
from pathlib import Path
import numpy as np, pandas as pd, pingouin as pg
from scipy import stats

MF = Path(data_dir + r"\Revision_2026\Manuscript_Figures"); sys.path.insert(0, str(MF / "code"))
import tep_io as io  # noqa: E402
RV = MF / "Review_2026-11-10"; CK = RV / "_check"; DATA = MF / "data"
LEV, EYES, S17 = io.LEV, io.EYES, io.MAIN17
WIN_TEP = {"P30": (28, 48), "N45": (41, 61), "P60": (48, 68), "N100": (120, 140)}
WIN_ERM = {"P30": (27, 47), "N45": (39, 59), "P60": (50, 70), "N100": (120, 140)}

_ERM = {}
def erm(s, e, stim="Real"):
    k = (s, e, stim)
    if k not in _ERM:
        import scipy.io as sio
        _ERM[k] = sio.loadmat(io.HHSA / f"erm_spTMS_{s:03d}_{stim}_{e}.mat", squeeze_me=True)
    return _ERM[k]
def erm_wave(s, e, lv, stim="Real"):
    M = erm(s, e, stim); return np.round(M["time"]).astype(int), np.asarray(M["erm_roi"][LEV.index(lv)], float)

_W = pd.read_csv(DATA / "TEP_roi_waveforms.csv"); _W = _W[_W.file == "TEP"]; _tc = [c for c in _W.columns if c.startswith("t_")]
T_TEP = np.arange(-2000, -2000 + len(_tc))
TEPW = {(int(r.id[6:9]), r.id.split("_")[2], r.id[-2:], int(r.intensity)): np.asarray(r[_tc], float) for _, r in _W.iterrows()}
def tep_wave(s, e, lv, stim="Real"): return T_TEP, TEPW[(s, stim, e, lv)]

def wmean(t, y, a0, a1): return y[(t >= a0) & (t <= a1)].mean()
def winmean21(t, y, a): i0 = np.where(t >= a)[0][0]; return y[i0:i0 + 21].mean()   # tep_stats.py sliding window (21 samples)

def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p); q = np.empty(m)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]; return np.minimum(q, 1)
def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p); adj = np.empty(m)
    adj[o] = np.minimum(np.maximum.accumulate((m - np.arange(m)) * p[o]), 1); return adj

def anova2(df, dv="y", within=("eye", "intensity"), subject="sub"):
    """2-way rmANOVA, pingouin, GG-corrected p; df reported as GG-corrected (df x epsilon)."""
    a = pg.rm_anova(data=df, dv=dv, within=list(within), subject=subject, correction=True, effsize="np2").set_index("Source"); out = []
    for src in a.index:
        r = a.loc[src]; eps = r["eps"] if "eps" in r and pd.notna(r.get("eps")) else 1.0
        pgg = r["p_GG_corr"] if "p_GG_corr" in r and pd.notna(r.get("p_GG_corr")) else r["p_unc"]
        out.append(dict(effect=src.replace(" * ", " x "), F=r["F"], df1=r["ddof1"] * eps, df2=r["ddof2"] * eps, p=pgg, eta_p2=r["np2"], eps=eps))
    return pd.DataFrame(out)

def paired(a, b):
    d = np.asarray(a, float) - np.asarray(b, float); r = stats.ttest_rel(a, b)
    return dict(n=len(d), mean_diff=d.mean(), SD_diff=d.std(ddof=1), dz=d.mean() / d.std(ddof=1), t=r.statistic, df=len(d) - 1, p=r.pvalue)

def fp(p):
    """APA p: no leading zero, three decimals, p < .001."""
    if p < .001: return "p < .001"
    s = f"{p:.3f}"; return "p = " + (s[1:] if s.startswith("0") else s)
def f2(x): return f"{x:.2f}".replace("-", "\u2212")
def lz(x, d=2):   # no leading zero (correlations, eta)
    s = f"{x:.{d}f}"; s = s.replace("0.", ".", 1) if s.startswith(("0.", "-0.")) else s; return s.replace("-", "\u2212")
def fF(r): return f"F({r['df1']:.2f}, {r['df2']:.2f}) = {r['F']:.2f}"

def md5(p):
    import hashlib; return hashlib.md5(Path(p).read_bytes()).hexdigest()
