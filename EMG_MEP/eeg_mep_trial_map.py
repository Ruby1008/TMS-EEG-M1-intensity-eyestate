"""EEG epoch <-> MEP trial mapping for the single-trial mixed model (master plan F), 17 participants, real stimulation.
Two candidate keys, both checked against the intensity of each MEP trial (VETA TMS_Condition 1/2/3 = 80/100/120% RMT):
  A) stimulation order: pulse ordinal (EEG.event.urevent, the order of the detected pulses) shifted by the stimuli that produced no
     pulse (data/trials_without_pulse.csv, stimOrderIndex) -> MEP trial_number (stimulation order; HANDOFF 一-2)
  B) trial2023 (EEG.event.trial2023; master plan F wording) -> MEP trial_number
Output: data/EEG_MEP_trial_map.csv (id, epoch, intensity, urevent, trial2023, stimorder_A, mep_trial_A/B, intensities, MEP amplitude, accept),
        data/EEG_MEP_trial_map_check.csv (per recording: agreement of intensities and number of usable trials for A and B)
"""
import numpy as np, pandas as pd
import tep_io as io

MEP = {e: pd.read_csv(io.ROOT / "MEP_source" / f"MEP_All23_{e}_20240416.csv", encoding="utf-8-sig") for e in io.EYES}
NOP = pd.read_csv(io.DATA / "trials_without_pulse.csv")
rows, chk = [], []
for s in io.MAIN17:
    for e in io.EYES:
        id_ = f"spTMS_{s:03d}_Real_{e}"
        _, _, _, inten, ev = io.load_set(io.FINAL, id_ + "_F_ERM")
        m = MEP[e][MEP[e]["sub"] == s].set_index("trial_number")
        gaps = sorted(NOP[NOP.id == id_].stimOrderIndex.tolist())
        def order(p):                      # p-th delivered pulse -> stimulation order, skipping the stimuli without a pulse
            o = p
            for g in gaps:
                if g <= o: o += 1
            return o
        for k, ev_ in enumerate(ev):
            ep = int(ev_.epoch); ue = int(ev_.urevent); t23 = int(ev_.trial2023) if np.size(ev_.trial2023) and not np.isnan(float(ev_.trial2023)) else -1
            oa = order(ue)
            r = dict(id=id_, sub=s, eye=e, epoch=ep, intensity=int(ev_.intensity), urevent=ue, trial2023=t23, stimorder_A=oa)
            for key, tr in (("A", oa), ("B", t23)):
                if tr in m.index:
                    mm = m.loc[tr]; r[f"mep_int_{key}"] = {1: 80, 2: 100, 3: 120}.get(int(mm.TMS_Condition), np.nan)
                    r[f"mep_amp_{key}"] = mm.ch4_MEP_amplitude; r[f"mep_accept_{key}"] = int(mm.trial_accept)
            rows.append(r)
        d = pd.DataFrame([r for r in rows if r["id"] == id_])
        c = dict(id=id_, n_epochs=len(d))
        for key in ("A", "B"):
            ok = d.get(f"mep_int_{key}")
            c[f"agree_{key}"] = float((d.intensity == d[f"mep_int_{key}"]).mean()) if ok is not None else np.nan
            c[f"usable_{key}"] = int(((d.intensity == d[f"mep_int_{key}"]) & (d[f"mep_accept_{key}"] == 1)).sum()) if ok is not None else 0
        chk.append(c)
pd.DataFrame(rows).to_csv(io.DATA / "EEG_MEP_trial_map.csv", index=False)
C = pd.DataFrame(chk); C.to_csv(io.DATA / "EEG_MEP_trial_map_check.csv", index=False)
pd.set_option("display.width", 200); print(C.round(3).to_string())
print("median agreement A", C.agree_A.median(), "min", C.agree_A.min(), "| B", C.agree_B.median(), "min", C.agree_B.min())
