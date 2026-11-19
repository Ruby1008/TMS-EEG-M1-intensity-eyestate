"""Per-condition all-channel averages for the supplementary TEP measures (master plan C; S-E band-pass control, decisions 2026-10-16).
All available real-stimulation recordings (stage-D exclusions removed). For each recording x intensity:
  TEP = _F_TEP average (0.1-30 Hz); RAW = _F_ERM average (0.1-80 Hz, no low-pass);
  BP  = RAW average band-passed 4.2-36.4 Hz (scipy butter(4, [4.2, 36.4], 'band'), filtfilt, zero phase; as the S-C control in
        erm_prestim_variants.py; filtering the trial average equals averaging the filtered trials because the filter is linear).
Saved: data/ALLCH_avg.npz (keys '<id>|<kind>' -> chan x time x 3 intensities, float32), labels, time (ms, -2000..1999).
"""
import numpy as np
from scipy import signal
import tep_io as io

b, a = signal.butter(4, [4.2, 36.4], btype="band", fs=1000)
out = {}; labs = None; t = None
ids = [f"spTMS_{s:03d}_Real_{e}" for e in io.EYES for s in io.AVAIL[e]]
for id_ in ids:
    for kind, suf in (("TEP", "_F_TEP"), ("RAW", "_F_ERM")):
        x, l, tt, inten, _ = io.load_set(io.FINAL, id_ + suf)
        if labs is None: labs, t = l, tt
        assert l == labs
        avg = np.stack([x[:, :, inten == lv].mean(2) for lv in io.LEV], axis=2).astype(np.float32)
        out[f"{id_}|{kind}"] = avg
        if kind == "RAW": out[f"{id_}|BP"] = signal.filtfilt(b, a, avg.astype(float), axis=1).astype(np.float32)
    print(id_)
np.savez_compressed(io.DATA / "ALLCH_avg.npz", labels=np.array(labs), time=t, **out)
print(len(ids), "recordings")
