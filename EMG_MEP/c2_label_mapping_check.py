"""Part 4 C2 (user 2026-10-23): label and mapping check for all 46 Real MEP sessions (method of qc/a1_022_mapping_check.py).
Per session:
 (a) raw trigger code (EEG Preproc/Intermediate/<id>_e.set, event types 1/2/3 = 80/100/120% RMT) vs VETA label (TMS_Condition) vs EEG label
     (data/EEG_MEP_trial_map.csv 'intensity' = label of the final EEG epoch), under the sequential mapping
     (_e epoch i = i-th VETA trial that is not in no_pulse.csv) - only valid when the counts agree;
 (b) key B (trial2023) vs sequential mapping for the final EEG epochs (linked through urevent): number of epochs with a different VETA trial,
     code agreement under both mappings, and final epochs mapped by key B to VETA trials without a pulse (no_pulse.csv).
Own code of an epoch = last code event 3-100 ms before the pulse (as pair_tms_intensity, maxLag 100; code events 0-1 ms around the pulse are
pulse artifacts on the trigger line, e.g. 022 EO epoch 44: '2' at -23 ms, '3' at +1 ms); code values per session from datalist evCodes
(123: 1/2/3, 246: 2/4/6 = 80/100/120% RMT). Sessions without an _e file (EEG excluded before processing) are listed with the reason.
Outputs: qc/label_mapping_check_all_sessions.csv (one row per session), qc/label_mapping_check_trials.csv (every disagreeing trial),
         qc/EEG_MEP_trial_map_C2fix.csv (key B replaced by the sequential mapping where key B is shifted; original in trial2023_keyB_orig)."""
data_dir = r"path/to/data"   # release: set this path
from pathlib import Path
import numpy as np, pandas as pd, scipy.io as sio
REV = Path(data_dir + r"\Revision_2026"); MF = REV / "Manuscript_Figures"; R = MF / "MEP_reprocess_2026-10-22"; Q = R / "qc"
I = pd.read_csv(R / "input_files.csv", dtype=str); NP = pd.read_csv(R / "no_pulse.csv", dtype=str)
meta = pd.read_csv(R / "extract" / "mep_trials_meta.csv", dtype={"subject": str}); meta["veta_label"] = meta.TMS_Condition.map({1: 80, 2: 100, 3: 120})
M = pd.read_csv(MF / "data" / "EEG_MEP_trial_map.csv"); M2 = M.copy(); M2["trial2023_keyB_orig"] = M2.trial2023
DL = pd.read_csv(REV / "Preproc" / "datalist_Real.csv", encoding="utf-8-sig").set_index("id")
CODE = {"1": 80, "2": 100, "3": 120}
rows, trows = [], []
for s, e in I[["subject", "eye"]].drop_duplicates().itertuples(index=False):
    sid = f"spTMS_{s}_Real_{e}"; ef = REV / "Preproc" / "Intermediate" / f"{sid}_e.set"
    v = meta[(meta.subject == s) & (meta.eye == e)].sort_values("trial_number").set_index("trial_number")
    nop = sorted(int(x) for x in NP[(NP.subject == s) & (NP.eye == e)].trial_number)
    pulsed = [n for n in v.index if n not in nop]
    r = dict(session=sid, veta_trials=len(v), no_pulse=len(nop), veta_pulsed=len(pulsed), evCodes=str(DL.loc[sid, "evCodes"]) if sid in DL.index else "", eeg_label_source=DL.loc[sid, "labelSource"] if sid in DL.index else "")
    if not ef.exists():
        r.update(status="no _e file (EEG excluded before processing)"); rows.append(r); continue
    E = sio.loadmat(str(ef), squeeze_me=True, struct_as_record=False); ntr = int(E["trials"])
    ev = str(DL.loc[sid, "evCodes"]); CM = {"123": {"1": 80, "2": 100, "3": 120}, "246": {"2": 80, "4": 100, "6": 120}}[ev]
    tms = {x.epoch: x.latency for x in E["event"] if x.type == "TMS"}; ure = {x.epoch: x.urevent for x in E["event"] if x.type == "TMS"}
    code = {}
    for x in E["event"]:                      # own code = last code event 3-100 ms before the pulse (maxLag 100 of pair_tms_intensity; code events at 0-1 ms are artifacts of the pulse on the trigger line)
        if x.type in CM and 3 <= tms[x.epoch] - x.latency <= 100 and (x.epoch not in code or x.latency > code[x.epoch][1]): code[x.epoch] = (CM[x.type], x.latency)
    code = {k: c for k, (c, _) in code.items()}
    EE = pd.DataFrame(dict(e_epoch=range(1, ntr + 1))); EE["code"] = EE.e_epoch.map(code); EE["urevent"] = EE.e_epoch.map(ure)
    r.update(e_epochs=ntr, epochs_without_code=int(EE.code.isna().sum()), count_match=ntr == len(pulsed))
    if ntr == len(pulsed):
        EE["veta_seq"] = pulsed; EE["veta_label_seq"] = EE.veta_seq.map(v.veta_label)
        ok = EE.code.notna(); r.update(code_vs_veta_seq=f"{int((EE.code[ok] == EE.veta_label_seq[ok]).sum())}/{int(ok.sum())}")
        for x in EE[ok & (EE.code != EE.veta_label_seq)].itertuples():
            trows.append(dict(session=sid, check="(a) code vs VETA label (sequential mapping)", e_epoch=x.e_epoch, veta_trial=x.veta_seq, code=x.code, veta_label=x.veta_label_seq, eeg_label=np.nan))
    else:
        r.update(code_vs_veta_seq="counts differ: sequential mapping not defined")
    m = M[M.id == sid]
    if len(m):
        u = EE.set_index("urevent"); mm = m.assign(code=m.urevent.map(u.code), e_epoch=m.urevent.map(u.e_epoch))
        mm["veta_keyB_label"] = mm.trial2023.map(v.veta_label)
        okc = mm.code.notna()
        r.update(final_epochs=len(m), eeg_label_vs_code=f"{int((mm.intensity[okc] == mm.code[okc]).sum())}/{int(okc.sum())}",
                 code_vs_veta_keyB=f"{int((mm.code[okc] == mm.veta_keyB_label[okc]).sum())}/{int(okc.sum())}", keyB_to_no_pulse=int(mm.trial2023.isin(nop).sum()))
        for x in mm[okc & (mm.intensity != mm.code)].itertuples():
            trows.append(dict(session=sid, check="(a) EEG label vs code", e_epoch=x.e_epoch, final_epoch=x.epoch, veta_trial=np.nan, code=x.code, veta_label=np.nan, eeg_label=x.intensity))
        if ntr == len(pulsed):
            mm["veta_seq"] = mm.urevent.map(u.veta_seq); mm["veta_seq_label"] = mm.veta_seq.map(v.veta_label); d = mm[mm.veta_seq != mm.trial2023]
            r.update(keyB_differs_from_seq=len(d), code_vs_veta_seq_final=f"{int((mm.code[okc] == mm.veta_seq_label[okc]).sum())}/{int(okc.sum())}")
            if len(d):
                sh = (d.veta_seq - d.trial2023).value_counts().to_dict(); r.update(keyB_shift=str({int(k): int(n) for k, n in sh.items()}), keyB_shift_range=f"{int(d.trial2023.min())}-{int(d.trial2023.max())}")
                for x in d.itertuples():
                    trows.append(dict(session=sid, check="(b) key B differs from sequential", e_epoch=x.e_epoch, final_epoch=x.epoch, veta_trial=x.veta_seq, keyB_trial=x.trial2023,
                                      code=x.code, veta_label=x.veta_seq_label, veta_label_keyB=x.veta_keyB_label, eeg_label=x.intensity))
                # correct only when the sequential mapping agrees better with the codes than key B
                if (mm.code[okc] == mm.veta_seq_label[okc]).sum() > (mm.code[okc] == mm.veta_keyB_label[okc]).sum():
                    M2.loc[M2.id == sid, "trial2023"] = mm.veta_seq.values.astype(int); r.update(mapping_corrected=True)
                else: r.update(mapping_corrected=False)
    r.update(status="checked"); rows.append(r)
S = pd.DataFrame(rows); S.to_csv(Q / "label_mapping_check_all_sessions.csv", index=False, encoding="utf-8-sig")
T = pd.DataFrame(trows); T.to_csv(Q / "label_mapping_check_trials.csv", index=False, encoding="utf-8-sig")
M2.to_csv(Q / "EEG_MEP_trial_map_C2fix.csv", index=False)
pd.set_option("display.width", 260); pd.set_option("display.max_columns", 30); print(S.to_string(index=False)); print(T.groupby(["session", "check"]).size().to_string())
