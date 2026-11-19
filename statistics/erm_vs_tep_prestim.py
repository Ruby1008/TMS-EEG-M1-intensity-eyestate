"""Section 3-1 (Supplementary): ERM (CEEMDAN IMF4-7, from HHSA/erm_<id>.mat) vs the conventional TEP
(_F_TEP, Stats_R1-1/prestim_phaseA_def.csv) on the same datasets.
  (a) overlay of the grand-average waveforms (ROI, per intensity, Real EO/EC)
  (b) pre-stimulus check on the ERM with the R1-1 method (Phase A definition: ref -300..-100, shift
      -20..-5 ms): shift == 0 (t, Wilcoxon, Holm), TOST +-0.5/+-1 uV (Holm), rmcorr(shift, 15-50 ms)
  (c) paired comparison ERM vs TEP of |shift| per dataset x intensity (is the ERM pre-stimulus offset larger?)
Stop rule (user 2026-10-12): if the ERM pre-stimulus offset is clearly larger than the TEP one, report before
running anything else.
Usage: python erm_vs_tep_prestim.py [--ids id1,id2]   (default: every erm_*.mat present)
"""
import sys, base64, io
from pathlib import Path
import numpy as np, pandas as pd, scipy.io as sio
from scipy import stats
import pingouin as pg
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import os, sys as _s; _s.path.insert(0, str(Path(__file__).resolve().parent)); from tep_io import MAIN17
SUB17 = os.environ.get("SUBSET") == "main17"; SFX = "_main17" if SUB17 else ""   # 2026-10-16: 17-participant main version

ROOT = Path(__file__).resolve().parent.parent; HH = ROOT / "HHSA"; OUT = ROOT / "data"; OUT.mkdir(exist_ok=True)
PRE = ROOT.parent / "Preproc"
BOUNDS = [0.5, 1.0]
ids = None
if "--ids" in sys.argv: ids = sys.argv[sys.argv.index("--ids") + 1].split(",")
files = sorted(HH.glob("erm_spTMS_*.mat")); files = [f for f in files if ids is None or f.stem[4:] in ids]
assert files, "no erm_*.mat files"

# --- ERM per dataset
rows, wave = [], {}
for f in files:
    R = sio.loadmat(f, squeeze_me=True)
    id_ = str(R["id"]); t = R["time"]; erm = R["erm_roi"]; wave[id_] = (t, erm)
    for k, lv in enumerate([80, 100, 120]):
        w = erm[k] - erm[k][(t >= -300) & (t <= -100)].mean()
        s = (t >= -60) & (t <= -5); pf = np.polyfit(t[s], w[s], 1)
        rows.append({"id": id_, "intensity": lv, "nTrials": int(R["nTrials"][k]), "shift": w[(t >= -20) & (t <= -5)].mean(),
                     "slope": pf[0] * 10, "post": w[(t >= 15) & (t <= 50)].mean()})
E = pd.DataFrame(rows)
# variant (a) (user 2026-10-12): ERM baseline-corrected to -100..-10 ms, then the same shift definition
# (shift = mean(-20..-5) after ref -300..-100) -> identical to the original by construction; what changes with a
# baseline correction is the ABSOLUTE offset, so we also report the raw offset mean(-20..-5) without the -300..-100 ref
rows_raw = []
for id_, (t, erm) in wave.items():
    for k, lv in enumerate([80, 100, 120]):
        w = erm[k]; rows_raw.append({"id": id_, "intensity": lv, "shift_raw": w[(t >= -20) & (t <= -5)].mean(),
                                     "shift_bl100": (w - w[(t >= -100) & (t <= -10)].mean())[(t >= -20) & (t <= -5)].mean()})
E = E.merge(pd.DataFrame(rows_raw), on=["id", "intensity"]); E["sub"] = E.id.str.extract(r"spTMS_(\d+)_").astype(int); E["stim"] = E.id.str.split("_").str[2]; E["eyes"] = E.id.str[-2:]
# --- TEP (same definition, _F_TEP) and stage-D set
T = pd.read_csv(ROOT.parent / "Preproc" / "Stats_R1-1" / "prestim_phaseA_def.csv"); T = T[T.file == "TEP"]
Q = pd.read_csv(PRE / "QC" / "QC_summary_long.csv", encoding="utf-8-sig")[["id", "intensity", "stageD"]]
M = E.merge(T[["id", "intensity", "shift_m20_m5", "slope_m60_m5_uVper10ms", "post_15_50"]], on=["id", "intensity"], how="left") \
     .rename(columns={"shift_m20_m5": "shift_TEP", "slope_m60_m5_uVper10ms": "slope_TEP", "post_15_50": "post_TEP"}).merge(Q, on=["id", "intensity"], how="left")
M.to_csv(OUT / f"S1_ERM_vs_TEP_prestim_per_dataset{SFX}.csv", index=False)
main = M[M.stageD != "exclude"]
if SUB17: main = main[main["sub"].isin(MAIN17)]
print(f"datasets: {M.id.nunique()} (Real {M[M.stim=='Real'].id.nunique()}, Sham {M[M.stim=='Sham'].id.nunique()}); main-set rows {len(main)}")

# --- (c) ERM vs TEP |shift|, paired
d = main.dropna(subset=["shift_TEP"])
w = stats.wilcoxon(d["shift"].abs(), d.shift_TEP.abs())
print(f"|shift| ERM median {d['shift'].abs().median():.3f} vs TEP {d.shift_TEP.abs().median():.3f} uV; "
      f"ERM larger in {(d['shift'].abs() > d.shift_TEP.abs()).mean()*100:.0f}% of {len(d)} cells; Wilcoxon p = {w.pvalue:.3g}")
print("raw offset mean(-20..-5) without ref: ERM |raw| median %.2f uV; after -100..-10 baseline |shift| median %.2f uV" % (d.shift_raw.abs().median(), d.shift_bl100.abs().median()))
print(f"shift SD across cells: ERM {d['shift'].std():.3f}, TEP {d.shift_TEP.std():.3f}; corr(ERM, TEP shift) r = {np.corrcoef(d['shift'], d.shift_TEP)[0,1]:.2f}")

# --- (b) R1-1 tests on the ERM shift, main set
def holm(p): return pg.multicomp(np.asarray(p, float), method="holm")[1]
def tost_p(x, b):
    n = len(x); m = x.mean(); se = x.std(ddof=1) / np.sqrt(n)
    return max(stats.t.sf((m + b) / se, n - 1), stats.t.cdf((m - b) / se, n - 1))
cells = [(s, e, i) for s in ("Real", "Sham") for e in ("EO", "EC") for i in (80, 100, 120)]
res = []
for s, e, i in cells:
    x = main[(main.stim == s) & (main.eyes == e) & (main.intensity == i)]["shift"].to_numpy()
    if len(x) < 5: continue
    n = len(x); m = x.mean(); se = x.std(ddof=1) / np.sqrt(n); ci = stats.t.interval(0.90, n - 1, loc=m, scale=se)
    r = {"stim": s, "eyes": e, "intensity": i, "n": n, "mean": m, "ci90_lo": ci[0], "ci90_hi": ci[1],
         "p_t": stats.ttest_1samp(x, 0).pvalue, "p_wilcoxon": stats.wilcoxon(x).pvalue if np.any(x) else 1}
    for b in BOUNDS: r[f"p_tost_{b:g}"] = tost_p(x, b)
    res.append(r)
Rt = pd.DataFrame(res)
for c in ["p_t", "p_wilcoxon"] + [f"p_tost_{b:g}" for b in BOUNDS]: Rt[c + "_holm"] = holm(Rt[c])
Rt.to_csv(OUT / f"S1_ERM_prestim_shift_tests{SFX}.csv", index=False)
print(Rt.round(3)[["stim", "eyes", "intensity", "n", "mean", "ci90_lo", "ci90_hi", "p_t_holm", "p_tost_0.5_holm", "p_tost_1_holm"]].to_string())
rm = []
for s in ("Real", "Sham"):
    dd = main[main.stim == s]
    if dd["sub"].nunique() < 5: continue
    r = pg.rm_corr(data=dd, x="shift", y="post", subject="sub"); ci = r["CI95" if "CI95" in r else "CI95%"].iloc[0]
    rm.append({"stim": s, "nSubjects": dd["sub"].nunique(), "nObs": len(dd), "r_rm": r["r"].iloc[0], "ci_lo": ci[0], "ci_hi": ci[1], "p": r["pval"].iloc[0]})
RM = pd.DataFrame(rm); RM.to_csv(OUT / f"S1_ERM_prestim_rmcorr{SFX}.csv", index=False); print(RM.round(3).to_string())

# --- (a) overlay figure: Real EO/EC x intensity, ERM vs TEP grand average (ROI)
def ga_tep(ids_):
    acc = {};
    return acc
TEPW = None; fw = OUT / "TEP_roi_waveforms.csv"
if fw.exists():
    TEPW = pd.read_csv(fw); tcols = [c for c in TEPW.columns if c.startswith("t_")]; tt = np.arange(-2000, -2000 + len(tcols))   # MATLAB mangles "t_-2000"; samples are -2000..1999 ms at 1 kHz
fig, axs = plt.subplots(2, 3, figsize=(13, 6), sharex=True, sharey=True)
cols = {80: "#f9954e", 100: "#f97213", 120: "#c55505"}
for r_, e in enumerate(("EO", "EC")):
    for c_, lv in enumerate((80, 100, 120)):
        ax = axs[r_, c_]; ids_ = [i for i in wave if i.endswith(f"Real_{e}") and i in set(main.id)]
        if not ids_: ax.set_title(f"{e} {lv}%: no data"); continue
        k = [80, 100, 120].index(lv); t = wave[ids_[0]][0]
        erm = np.mean([wave[i][1][k] for i in ids_], axis=0)
        ax.plot(t, erm, color=cols[lv], lw=2, label="ERM (IMF4-7)")
        if TEPW is not None:
            for fl, ls, lab in (("TEP", "-", "_F_TEP (30 Hz LP)"), ("ERM", ":", "_F_ERM (no LP)")):
                g = TEPW[(TEPW.file == fl) & (TEPW.intensity == lv) & TEPW.id.isin(ids_)]
                if len(g): ax.plot(tt, g[tcols].to_numpy().mean(0), color="k", lw=1.2, ls=ls, label=lab)
        ax.axvspan(-5, 10, color="k", alpha=.8); ax.axhline(0, color="#999", lw=.7); ax.set_xlim(-100, 400)
        ax.set_title(f"Real {e} {lv}% RMT"); ax.legend(fontsize=8, loc="upper right")
axs[1, 0].set_xlabel("Time (ms)"); axs[0, 0].set_ylabel("ROI (µV)"); axs[1, 0].set_ylabel("ROI (µV)")
fig.suptitle("Grand average, ROI FC3/FC1/C3/C1: ERM (CEEMDAN IMF4-7) vs conventional TEP", fontsize=11)
fig.tight_layout(); fig.savefig(ROOT / f"FigS_ERMvsTEP_overlay{SFX}.png", dpi=600); plt.close(fig)
print("wrote", OUT)
