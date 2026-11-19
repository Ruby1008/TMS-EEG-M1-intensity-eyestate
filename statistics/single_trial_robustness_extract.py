"""Single-trial variables for the robustness of P30 -> log MEP (user 2026-10-20, task 1). One row per EEG epoch of the 17 main participants x EO/EC (real).
  P30         : ERM (IMFs 4-7) ROI mean 27-47 ms (data/ERM_single_trial.csv; erm_single_trial_extract.m)
  P30_BP      : _F_ERM (0.1-80 Hz) ROI mean, band-passed 4.2-36.4 Hz (butter(4), filtfilt, per trial over the whole -2..2 s epoch; as the S-C/S-E control), 27-47 ms
  P30_indiv   : ERM ROI mean within +-10 ms of the participant's individual P30 latency (data/TEP_supp_indiv_latencies.csv, data = ERM; one latency per participant),
                single trials from HHSA\\timf_spm_<id> (IMFs 4-7, ROI FC3/FC1/C3/C1), same reading as erm_single_trial_windows.py
  resid_11_25 : maximum absolute ROI mean (_F_ERM, 0.1-80 Hz) 11-25 ms (the residual-artefact measure of the exclusion rule, per trial)
  trial_order : stimulation order within the block (data/EEG_MEP_trial_map.csv, stimorder_A)
  log_mu      : pre-stimulus mu power, data/Prestim_power_single_trial.csv (8-12 Hz, -1000..-10 ms, Welch 500-sample Hann, 50% overlap, log10, ROI)
  MEP fields  : trial2023 mapping (key B), as Figure 6
Output: data/Single_trial_robustness_input.csv"""
import numpy as np, pandas as pd, scipy.io as sio
from scipy import signal
import tep_io as io

b, a = signal.butter(4, [4.2, 36.4], btype="band", fs=1000)
ST = pd.read_csv(io.DATA / "ERM_single_trial.csv"); MP = pd.read_csv(io.DATA / "EEG_MEP_trial_map.csv"); PW = pd.read_csv(io.DATA / "Prestim_power_single_trial.csv")
LAT = pd.read_csv(io.DATA / "TEP_supp_indiv_latencies.csv"); LAT = LAT[(LAT.data == "ERM") & (LAT.component == "P30")].set_index("sub").latency_ms
rows = []
for s in io.MAIN17:
    lat = int(LAT.loc[s])
    for e in io.EYES:
        id_ = f"spTMS_{s:03d}_Real_{e}"
        x, labs, t, inten, _ = io.load_set(io.FINAL, id_ + "_F_ERM"); t = np.round(t).astype(int); roi = [labs.index(c) for c in io.ROI]
        R = x[roi].astype(float).mean(0)                                                   # time x trial
        Rbp = signal.filtfilt(b, a, R, axis=0)
        w30 = (t >= 27) & (t <= 47); wres = (t >= 11) & (t <= 25)
        p30_bp = Rbp[w30].mean(0); res = np.abs(R[wres]).max(0)
        S = sio.loadmat(io.HHSA / f"timf_spm_{id_}.mat", squeeze_me=True, struct_as_record=False)["D"]
        nch, nimf, nt, ntr = (int(v) for v in S.data.dim); lab2 = [c.label for c in np.atleast_1d(S.channels)]; r2 = [lab2.index(c) for c in io.ROI]
        tt = np.round((S.timeOnset + np.arange(int(S.Nsamples)) / float(S.Fsample)) * 1000).astype(int); w0 = np.where(tt == lat - 10)[0][0]; w1 = np.where(tt == lat + 10)[0][0]
        inten2 = np.array([float(tr.label) for tr in np.atleast_1d(S.trials)]); assert ntr == len(inten) and np.array_equal(inten2, inten)
        p30_ind = np.empty(ntr); blk = nch * nimf
        with open(io.HHSA / f"timf_spm_{id_}.dat", "rb") as fh:
            for k in range(ntr):
                fh.seek(4 * (blk * nt * k + blk * w0)); A = np.fromfile(fh, dtype="<f4", count=blk * (w1 - w0 + 1)).reshape((nch, nimf, w1 - w0 + 1), order="F")
                p30_ind[k] = A[r2][:, 3:7].mean(0).sum(0).mean()
        for k in range(ntr):
            rows.append(dict(id=id_, epoch=k + 1, intensity=int(inten[k]), P30_BP=p30_bp[k], P30_indiv=p30_ind[k], indiv_latency=lat, resid_11_25=res[k]))
        print(id_, ntr, flush=True)
X = pd.DataFrame(rows)
X = X.merge(ST[["id", "epoch", "intensity", "P30"]], on=["id", "epoch"], suffixes=("", "_st")); assert (X.intensity == X.intensity_st).all()
X = X.merge(MP[["id", "epoch", "sub", "eye", "stimorder_A", "mep_int_B", "mep_amp_B", "mep_accept_B"]], on=["id", "epoch"])
X = X.merge(PW[["id", "epoch", "log_mu"]], on=["id", "epoch"]); X = X.rename(columns={"stimorder_A": "trial_order"}).drop(columns="intensity_st")
assert len(X) == len(ST)
X.to_csv(io.DATA / "Single_trial_robustness_input.csv", index=False); print(len(X), "trials; corr(P30, P30_BP) =", round(X[["P30", "P30_BP"]].corr().iloc[0, 1], 3), "; corr(P30, P30_indiv) =", round(X[["P30", "P30_indiv"]].corr().iloc[0, 1], 3))
