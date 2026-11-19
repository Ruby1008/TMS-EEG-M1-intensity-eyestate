"""Shared readers for the TEP analyses (2026-10-16): EEGLAB .set (two-file) and erm_<id>.mat; the 17 main participants."""
from pathlib import Path
import numpy as np, pandas as pd, scipy.io as sio

ROOT = Path(__file__).resolve().parents[1]; REV = ROOT.parent; DATA = ROOT / "data"
FINAL = REV / "Preproc" / "Final"; HHSA = ROOT / "HHSA"
ROI = ["FC3", "FC1", "C3", "C1"]; LEV = [80, 100, 120]; EYES = ["EO", "EC"]
Q = pd.read_csv(REV / "Preproc" / "QC" / "QC_summary_wide.csv", encoding="utf-8-sig")
_m = Q[(Q.status == "processed") & (Q.stim == "Real") & (Q.stageD != "exclude")]
AVAIL = {e: sorted(_m[_m.eyes == e]["sub"]) for e in EYES}                 # robustness: all available per condition
MAIN17 = sorted(set(AVAIL["EO"]) & set(AVAIL["EC"])); assert len(MAIN17) == 17


def load_set(folder, name):
    """EEGLAB .set -> data (chan x time x trial, float32), labels, times (ms), intensity per trial, event list"""
    S = sio.loadmat(Path(folder) / f"{name}.set", squeeze_me=True, struct_as_record=False)
    E = S.get("EEG"); get = (lambda k: getattr(E, k)) if E is not None else (lambda k: S[k])
    nb, npn, ntr = int(get("nbchan")), int(get("pnts")), int(get("trials"))
    dat = get("data")
    x = np.fromfile(Path(folder) / dat, dtype="<f4").reshape((nb, npn, ntr), order="F") if isinstance(dat, str) else np.asarray(dat, float)
    labs = [c.labels for c in np.atleast_1d(get("chanlocs"))]
    ev = np.atleast_1d(get("event")); inten = np.full(ntr, np.nan)
    for e in ev: inten[int(e.epoch) - 1] = float(e.type)
    return x, labs, np.asarray(get("times"), float), inten, ev


def load_erm(sid, eye, folder=HHSA):
    return sio.loadmat(Path(folder) / f"erm_spTMS_{sid:03d}_Real_{eye}.mat", squeeze_me=True)
