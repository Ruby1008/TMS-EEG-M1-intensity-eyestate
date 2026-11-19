"""N5 extraction (user 2026-10-28). Real stimulation, 17 participants, every epoch (same epochs as data/ERM_single_trial.csv and
data/Prestim_power_single_trial.csv; checked by count and intensity per recording).
P30 variants per trial (ROI FC3/FC1/C3/C1):
  P30_main   = data/ERM_single_trial.csv P30 (ERM IMF4-7, 27-47 ms)
  P30_IMF45  = ROI_IMF4_P30 + ROI_IMF5_P30 of data/ERM_single_trial_windows.csv (ERM is the sum of IMFs 4-7; the sum of the four IMF columns equals P30_main)
  P30_IMF6   = ROI_IMF6_P30 of data/ERM_single_trial_windows.csv
  P30_convTEP = mean of the conventional TEP (Preproc/Final/<id>_F_TEP.set, 0.1-30 Hz) 28-48 ms, ROI
Predictor (e): IMF6 pre-stimulus amplitude envelope = |Hilbert transform| of the single-trial IMF6 (HHSA/timf_spm_<id>, first-layer CEEMDAN) computed
  per ROI channel on -2000..-1 ms, mean over -500..-100 ms, mean over the four ROI channels (µV); log10 of it also stored.
Output: _check/N5_single_trial.csv"""
import numpy as np, pandas as pd, scipy.io as sio
from scipy.signal import hilbert
import r1110_lib as L
io = L.io
ST = pd.read_csv(L.DATA / "ERM_single_trial.csv"); SW = pd.read_csv(L.DATA / "ERM_single_trial_windows.csv"); PW = pd.read_csv(L.DATA / "Prestim_power_single_trial.csv")
D = ST.merge(SW[["id", "epoch", "ROI_IMF4_P30", "ROI_IMF5_P30", "ROI_IMF6_P30", "ROI_IMF7_P30"]], on=["id", "epoch"]).merge(PW, on=["id", "epoch", "intensity"])
assert len(D) == len(ST) == len(PW)
D["P30_main"] = D.P30; D["P30_IMF45"] = D.ROI_IMF4_P30 + D.ROI_IMF5_P30; D["P30_IMF6"] = D.ROI_IMF6_P30
assert np.abs(D[["ROI_IMF4_P30", "ROI_IMF5_P30", "ROI_IMF6_P30", "ROI_IMF7_P30"]].sum(1) - D.P30).max() < 1e-4
conv, env = {}, {}
for s in L.S17:
    for e in L.EYES:
        id_ = f"spTMS_{s:03d}_Real_{e}"; st = ST[ST.id == id_].sort_values("epoch")
        x, labs, t, inten, _ = io.load_set(io.FINAL, id_ + "_F_TEP"); t = np.round(t).astype(int)
        assert len(inten) == len(st) and np.array_equal(inten.astype(int), st.intensity.values), id_
        ri = [labs.index(c) for c in io.ROI]; w = (t >= 28) & (t <= 48)
        v = x[ri][:, w, :].astype(float).mean((0, 1))
        for k in range(len(inten)): conv[(id_, k + 1)] = v[k]
        S = sio.loadmat(io.HHSA / f"timf_spm_{id_}.mat", squeeze_me=True, struct_as_record=False)["D"]
        nch, nimf, nt, ntr = (int(q) for q in S.data.dim); assert int(S.data.dtype) == 16 and ntr == len(inten)
        clabs = [c.label for c in np.atleast_1d(S.channels)]; tt = np.round((S.timeOnset + np.arange(int(S.Nsamples)) / float(S.Fsample)) * 1000).astype(int)
        cri = [clabs.index(c) for c in io.ROI]; i0, i1 = np.where(tt == -2000)[0][0], np.where(tt == -1)[0][0]; blk = nch * nimf
        sel = (tt[i0:i1 + 1] >= -500) & (tt[i0:i1 + 1] <= -100)
        with open(io.HHSA / f"timf_spm_{id_}.dat", "rb") as fh:
            for k in range(ntr):
                fh.seek(4 * (blk * nt * k + blk * i0))
                a = np.fromfile(fh, dtype="<f4", count=blk * (i1 - i0 + 1)).reshape((nch, nimf, i1 - i0 + 1), order="F")
                sig = a[cri, 5, :].astype(float)                                 # ROI channels, IMF6, -2000..-1 ms
                env[(id_, k + 1)] = np.abs(hilbert(sig, axis=1))[:, sel].mean()
        print(id_, ntr, flush=True)
D["P30_convTEP"] = [conv[(i, ep)] for i, ep in zip(D.id, D.epoch)]
D["env_IMF6"] = [env[(i, ep)] for i, ep in zip(D.id, D.epoch)]; D["log_env_IMF6"] = np.log10(D.env_IMF6)
D["sub"] = D.id.str.split("_").str[1].astype(int)
D[["id", "sub", "epoch", "intensity", "log_mu", "log_occ_alpha", "env_IMF6", "log_env_IMF6", "P30_main", "P30_IMF45", "P30_IMF6", "P30_convTEP"]].to_csv(L.CK / "N5_single_trial.csv", index=False)
print(len(D), "trials; corr P30_main with variants:", D[["P30_main", "P30_IMF45", "P30_IMF6", "P30_convTEP"]].corr().round(3).iloc[0].to_dict(),
      "| corr log_mu vs env:", round(np.corrcoef(D.log_mu, D.log_env_IMF6)[0, 1], 3))
