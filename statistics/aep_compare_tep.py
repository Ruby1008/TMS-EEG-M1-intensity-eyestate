"""AEP-ICA robustness (HANDOFF section 15; master plan J): conventional TEP before vs after AEP removal.
Before = Preproc\\Final\\<id>_F_TEP.set; after = Preproc_AEPr\\Final\\<id>_F_TEP_AEPr.set (aep_remove.m).
Participants: the 17 with complete real-stimulation EEG (Preproc\\QC\\QC_summary_wide.csv; master plan section 一).
ROI FC3/FC1/C3/C1, mean over trials of each intensity.
N100 = mean 120-140 ms (master plan: fixed window); 2 (eye) x 3 (intensity) rmANOVA, GG (pingouin), before and after.
P30/N45/P60 are added after the component windows are confirmed (stop point S2).
Outputs: data\\AEP_TEP_roi_waveforms.csv, data\\AEP_TEP_N100_values.csv, data\\AEP_TEP_N100_anova.csv,
         FigS_AEP_TEP_before_after.png / .tiff (600 dpi)
"""
from pathlib import Path
import numpy as np, pandas as pd, pingouin as pg, scipy.io as sio
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]; REV = ROOT.parent; D = ROOT / "data"
SRC = {"before": (REV / "Preproc" / "Final", "_F_TEP"), "after": (REV / "Preproc_AEPr" / "Final", "_F_TEP_AEPr")}
ROI = ["FC3", "FC1", "C3", "C1"]; LEV = [80, 100, 120]
Q = pd.read_csv(REV / "Preproc" / "QC" / "QC_summary_wide.csv", encoding="utf-8-sig")
m = Q[(Q.status == "processed") & (Q.stim == "Real") & (Q.stageD != "exclude")]
SUBS = sorted(set(m[m.eyes == "EO"]["sub"]) & set(m[m.eyes == "EC"]["sub"])); assert len(SUBS) == 17, SUBS

def load_set(folder, name):
    S = sio.loadmat(folder / f"{name}.set", squeeze_me=True, struct_as_record=False)
    E = S["EEG"] if "EEG" in S else None
    get = (lambda k: getattr(E, k)) if E is not None else (lambda k: S[k])
    nb, npn, ntr = int(get("nbchan")), int(get("pnts")), int(get("trials"))
    dat = get("data")
    if isinstance(dat, str):
        x = np.fromfile(folder / dat, dtype="<f4").reshape((nb, npn, ntr), order="F")
    else:
        x = np.asarray(dat, dtype=float)
    labs = [c.labels for c in np.atleast_1d(get("chanlocs"))]
    times = np.asarray(get("times"), float)
    ev = np.atleast_1d(get("event")); inten = np.full(ntr, np.nan)
    for e in ev: inten[int(e.epoch) - 1] = float(e.type)
    return x, labs, times, inten

rows, wav = [], []
for s in SUBS:
    for eye in ("EO", "EC"):
        id_ = f"spTMS_{s:03d}_Real_{eye}"
        for ver, (folder, suf) in SRC.items():
            x, labs, t, inten = load_set(folder, id_ + suf)
            ri = [labs.index(c) for c in ROI]; w = (t >= -100) & (t <= 400); nw = (t >= 120) & (t <= 140)
            for lv in LEV:
                y = x[ri][:, :, inten == lv].mean(2).mean(0)
                rows.append(dict(sub=s, eye=eye, intensity=lv, version=ver, N100=y[nw].mean(), n_trials=int((inten == lv).sum())))
                wav.append(pd.DataFrame(dict(sub=s, eye=eye, intensity=lv, version=ver, time_ms=t[w], uV=y[w])))
V = pd.DataFrame(rows); V.to_csv(D / "AEP_TEP_N100_values.csv", index=False)
W = pd.concat(wav); W.to_csv(D / "AEP_TEP_roi_waveforms.csv", index=False)
res = []
for ver in ("before", "after"):
    a = pg.rm_anova(data=V[V.version == ver], dv="N100", within=["eye", "intensity"], subject="sub", correction=True, effsize="np2")
    a.insert(0, "version", ver); a.insert(1, "n", len(SUBS)); res.append(a)
A = pd.concat(res); A.to_csv(D / "AEP_TEP_N100_anova.csv", index=False)
P = V.pivot_table(index=["sub", "eye", "intensity"], columns="version", values="N100").reset_index(); P["after_minus_before"] = P.after - P.before
desc = P.groupby(["eye", "intensity"]).agg(before_mean=("before", "mean"), after_mean=("after", "mean"), diff_mean=("after_minus_before", "mean"), diff_sd=("after_minus_before", "std")).reset_index()
desc.to_csv(D / "AEP_TEP_N100_descriptives.csv", index=False)

# figure: grand average before (solid) / after (dashed) / difference (grey), EO and EC rows x 3 intensities
plt.rcParams.update({"font.family": "Arial", "font.size": 9})
G = W.groupby(["eye", "intensity", "version", "time_ms"]).uV.mean().unstack("version").reset_index()
fig, axs = plt.subplots(2, 3, figsize=(11, 5.6), sharex=True, sharey=True)
for r, eye in enumerate(("EO", "EC")):
    for c, lv in enumerate(LEV):
        g = G[(G.eye == eye) & (G.intensity == lv)]; ax = axs[r, c]
        ax.axvspan(120, 140, color="0.9", lw=0); ax.axvspan(-5, 10, color="0.8", lw=0)
        ax.plot(g.time_ms, g.before, "k-", lw=1.2, label="Before AEP removal")
        ax.plot(g.time_ms, g.after, color="#c0392b", lw=1.2, ls="--", label="After AEP removal")
        ax.plot(g.time_ms, g.after - g.before, color="0.5", lw=1, label="After − before")
        ax.axhline(0, color="0.6", lw=.5); ax.axvline(0, color="0.6", lw=.5, ls=":")
        ax.set_title(f"Real {eye}, {lv}% RMT"); ax.set_xlim(-100, 400)
        if c == 0: ax.set_ylabel("ROI amplitude (µV)")
        if r == 1: ax.set_xlabel("Time (ms)")
axs[0, 0].legend(frameon=False, fontsize=7, loc="lower left")
fig.suptitle(f"Conventional TEP (0.1–30 Hz), ROI FC3/FC1/C3/C1, before and after AEP-ICA removal, n = {len(SUBS)}", fontsize=10)
fig.tight_layout()
fig.savefig(ROOT / "FigS_AEP_TEP_before_after.png", dpi=600); fig.savefig(ROOT / "FigS_AEP_TEP_before_after.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"}); plt.close(fig)
pd.set_option("display.width", 200)
print(A[["version", "n", "Source", "ddof1", "ddof2", "F", "p_unc", "p_GG_corr", "np2", "eps"]].round(4).to_string())
print(desc.round(3).to_string())
