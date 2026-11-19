"""Part 4 A1 (user 2026-10-22): 022 EO EEG epoch <-> VETA trial mapping check.
EEG: Preproc/Intermediate/spTMS_022_Real_EO_e.set (124 epochs, TMS artifact not yet removed). Per epoch: FC3 artifact = max |x - baseline| from -2 to 6 ms
     (baseline -100 to -10 ms); raw trigger code (event types 1/2/3 = 80/100/120% RMT); urevent of the TMS event (links to the final epochs and the map).
EMG: extract/mep_waveforms.csv (VETA export, -100..150 ms re trigger). Per VETA trial: artifact = max |y - median(-60..-30 ms)| from -30 to 2 ms; VETA label.
Mappings compared for the 116 final EEG epochs of 022 EO:
  key B (data/EEG_MEP_trial_map.csv, trial2023 = VETA trial_number), and
  sequential (EEG _e epoch i = i-th VETA trial with an EMG pulse, i.e. VETA 1-58 and 61-126; VETA 59 and 60 have no EMG artifact = no pulse).
Outputs (qc/): A1_022EO_trial_level.csv, A1_022EO_evidence.csv, EEG_MEP_trial_map_022fix.csv (map with the corrected trial2023 for 022 EO, original in
trial2023_keyB_orig), intensity_overrides_022fix.csv (spec rule: VETA label != EEG label under the corrected mapping -> EEG label)."""
data_dir = r"path/to/data"   # release: set this path
from pathlib import Path
import numpy as np, pandas as pd, scipy.io as sio
from scipy import stats
REV = Path(data_dir + r"\Revision_2026"); MF = REV / "Manuscript_Figures"; R = MF / "MEP_reprocess_2026-10-22"; Q = R / "qc"
E = sio.loadmat(str(REV / "Preproc" / "Intermediate" / "spTMS_022_Real_EO_e.set"), squeeze_me=True, struct_as_record=False)
nb, ntr, npt = int(E["nbchan"]), int(E["trials"]), int(E["pnts"])
X = np.fromfile(REV / "Preproc" / "Intermediate" / "spTMS_022_Real_EO_e.fdt", dtype="<f4").reshape((nb, npt, ntr), order="F")
c = [ch.labels for ch in E["chanlocs"]].index("FC3"); t = np.arange(npt) + E["xmin"] * 1000
code = {e.epoch: {"1": 80, "2": 100, "3": 120}[e.type] for e in E["event"] if e.type in ("1", "2", "3")}
ure = {e.epoch: e.urevent for e in E["event"] if e.type == "TMS"}
EE = pd.DataFrame([dict(e_epoch=i, urevent=ure[i], trigger_code=code[i],
                        eeg_fc3_artifact_uV=float(np.abs(X[c, (t >= -2) & (t <= 6), i - 1] - X[c, (t >= -100) & (t < -10), i - 1].mean()).max())) for i in range(1, ntr + 1)])
meta = pd.read_csv(R / "extract" / "mep_trials_meta.csv", dtype={"subject": str}); m = meta[(meta.subject == "022") & (meta.eye == "EO")].sort_values("trial_number")
W = pd.read_csv(R / "extract" / "mep_waveforms.csv"); W = W[W.row_id.isin(m.row_id)].set_index("row_id"); tw = np.array([int(s[1:]) for s in W.columns])
def emg_art(rid):
    y = W.loc[rid].values.astype(float); return float(np.abs(y[(tw >= -30) & (tw <= 2)] - np.median(y[(tw >= -60) & (tw <= -30)])).max())
V = m.assign(veta_label=m.TMS_Condition.map({1: 80, 2: 100, 3: 120}), emg_artifact_mV=[emg_art(r) for r in m.row_id]).set_index("trial_number")
NP = pd.read_csv(R / "no_pulse.csv", dtype=str); nop = sorted(int(x) for x in NP[(NP.subject == "022") & (NP.eye == "EO")].trial_number)
pulsed = [n for n in V.index if n not in nop]; assert len(pulsed) == ntr == 124 and nop == [59, 60]
EE["veta_sequential"] = pulsed
M = pd.read_csv(MF / "data" / "EEG_MEP_trial_map.csv"); sel = (M["sub"] == 22) & (M.eye == "EO"); m22 = M[sel].copy()
u = EE.set_index("urevent"); m22["e_epoch"] = m22.urevent.map(u.e_epoch); m22["trigger_code"] = m22.urevent.map(u.trigger_code); m22["eeg_fc3_artifact_uV"] = m22.urevent.map(u.eeg_fc3_artifact_uV)
m22["veta_keyB"] = m22.trial2023; m22["veta_sequential"] = m22.urevent.map(u.veta_sequential); m22["shift"] = m22.veta_sequential - m22.veta_keyB
for k in ("keyB", "sequential"):
    m22[f"veta_label_{k}"] = m22[f"veta_{k}"].map(V.veta_label); m22[f"emg_artifact_{k}"] = m22[f"veta_{k}"].map(V.emg_artifact_mV)
m22["emg_artifact_veta_59_60"] = np.nan
m22[["epoch", "intensity", "urevent", "e_epoch", "trigger_code", "eeg_fc3_artifact_uV", "veta_keyB", "veta_sequential", "shift", "veta_label_keyB", "veta_label_sequential",
     "emg_artifact_keyB", "emg_artifact_sequential"]].to_csv(Q / "A1_022EO_trial_level.csv", index=False)
ev = []
def add(item, keyB, seq, note=""): ev.append(dict(item=item, key_B=keyB, sequential=seq, note=note))
d = m22[m22["shift"] != 0]
add("final EEG epochs with a different VETA trial", len(d), len(d), f"key B {int(d.veta_keyB.min())}-{int(d.veta_keyB.max())}; shift -1 for {int((d['shift'] == -1).sum())}, +1 for {int((d['shift'] == 1).sum())}")
add("VETA label == EEG final label (116 epochs)", int((m22.veta_label_keyB == m22.intensity).sum()), int((m22.veta_label_sequential == m22.intensity).sum()))
add("VETA label == raw trigger code (116 epochs)", int((m22.veta_label_keyB == m22.trigger_code).sum()), int((m22.veta_label_sequential == m22.trigger_code).sum()))
add("VETA label == raw trigger code (27 differing epochs)", int((d.veta_label_keyB == d.trigger_code).sum()), int((d.veta_label_sequential == d.trigger_code).sum()))
add("EEG epochs mapped to VETA 59/60 (EMG without pulse artifact)", int(m22.veta_keyB.isin(nop).sum()), int(m22.veta_sequential.isin(nop).sum()),
    f"EMG artifact of VETA 59 / 60: {V.loc[59, 'emg_artifact_mV']:.4f} / {V.loc[60, 'emg_artifact_mV']:.4f} mV (median of all other trials {V.drop(nop).emg_artifact_mV.median():.4f})")
r_all = {k: stats.spearmanr(m22.eeg_fc3_artifact_uV, m22[f"emg_artifact_{k}"]) for k in ("keyB", "sequential")}
add("Spearman EEG FC3 artifact vs EMG artifact, 116 epochs", f"{r_all['keyB'].statistic:.3f} (p = {r_all['keyB'].pvalue:.4f})", f"{r_all['sequential'].statistic:.3f} (p = {r_all['sequential'].pvalue:.4f})")
r_d = {k: stats.spearmanr(d.eeg_fc3_artifact_uV, d[f"emg_artifact_{k}"]) for k in ("keyB", "sequential")}
add("Spearman EEG FC3 artifact vs EMG artifact, 27 differing epochs", f"{r_d['keyB'].statistic:.3f} (p = {r_d['keyB'].pvalue:.3f})", f"{r_d['sequential'].statistic:.3f} (p = {r_d['sequential'].pvalue:.3f})",
    "EMG artifact depends only weakly on intensity (see next rows), so this comparison has little power")
for k in (-2, -1, 0, 1, 2):
    vt = [pulsed[i + k] if 0 <= i + k < len(pulsed) else None for i in range(ntr)]; ea = np.array([V.loc[v, "emg_artifact_mV"] if v else np.nan for v in vt]); ok = ~np.isnan(ea)
    vl = np.array([V.loc[v, "veta_label"] if v else np.nan for v in vt])
    add(f"124 _e epochs, sequential mapping shifted by {k:+d}: Spearman EEG vs EMG artifact / label agreement", f"", f"{stats.spearmanr(EE.eeg_fc3_artifact_uV[ok], ea[ok]).statistic:.3f} / {np.nanmean(vl == EE.trigger_code):.3f}")
add("EEG FC3 artifact (median uV) by trigger code 80/100/120", "", " / ".join(f"{EE[EE.trigger_code == l].eeg_fc3_artifact_uV.median():.0f}" for l in (80, 100, 120)), f"Kruskal p = {stats.kruskal(*[g.eeg_fc3_artifact_uV for _, g in EE.groupby('trigger_code')]).pvalue:.1e}")
add("EMG artifact (median mV) by VETA label 80/100/120", "", " / ".join(f"{V[V.veta_label == l].emg_artifact_mV.median():.4f}" for l in (80, 100, 120)), f"Kruskal p = {stats.kruskal(*[g.emg_artifact_mV for _, g in V.drop(nop).groupby('veta_label')]).pvalue:.4f}")
EV = pd.DataFrame(ev); EV.to_csv(Q / "A1_022EO_evidence.csv", index=False, encoding="utf-8-sig")
# corrected map and overrides (spec section: VETA label differing from the EEG label -> EEG label)
M2 = M.copy(); M2["trial2023_keyB_orig"] = M2.trial2023; M2.loc[sel, "trial2023"] = m22.veta_sequential.values.astype(int)
M2.to_csv(Q / "EEG_MEP_trial_map_022fix.csv", index=False)
OV = m22.loc[m22.veta_label_sequential != m22.intensity, ["veta_sequential", "intensity"]].rename(columns={"veta_sequential": "trial_number"})
OV.insert(0, "eye", "EO"); OV.insert(0, "subject", "022"); OV.to_csv(Q / "intensity_overrides_022fix.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120); print(EV.to_string(index=False)); print(OV.to_string(index=False))
print(m22.loc[m22.veta_label_sequential != m22.intensity, ["epoch", "intensity", "trigger_code", "veta_sequential", "veta_label_sequential"]].to_string(index=False))
