"""Part 4 C2: raw .cnt trigger codes (qc/c2_cnt_events.csv, written by qc/c2_cnt_codes.m) vs VETA labels for the six sessions without an EEG _e file.
Codes 1/2/3 (= 80/100/120% RMT); a code event less than 1000 ms after the previous code event is treated as a duplicate of the same trial
(pulse artifact on the trigger line, as in the _e files where a second code follows the pulse by 0-1 ms). The code sequence is compared with
the 126 VETA labels in trial order, with offsets -3..+3. Every .cnt read (datalist file and the other recordings of 005, 006, 013) is compared with both VETA sessions of the participant.
Output: qc/label_mapping_check_cnt_sessions.csv"""
from pathlib import Path
import numpy as np, pandas as pd
R = Path(__file__).parents[1]; Q = R / "qc"
E = pd.read_csv(Q / "c2_cnt_events.csv", dtype={"type": str}, on_bad_lines="skip"); meta = pd.read_csv(R / "extract" / "mep_trials_meta.csv", dtype={"subject": str})
meta["lab"] = meta.TMS_Condition.map({1: 80, 2: 100, 3: 120}); CM = {"1": 80, "2": 100, "3": 120}
rows = []
for cnt_s in E.session.unique():
    e = E[(E.session == cnt_s) & E.type.isin(CM)].sort_values("latency_ms"); keep, last = [], -1e9
    for x in e.itertuples():
        if x.latency_ms - last >= 1000: keep.append(CM[x.type])
        last = x.latency_ms
    sub = cnt_s.split("_")[1]
    for eye in ("EO", "EC"):
        t = f"spTMS_{sub}_Real_{eye}"; v = meta[(meta.subject == sub) & (meta.eye == eye)].sort_values("trial_number").lab.tolist()
        best = (0, 0, 0)
        for k in range(-3, 4):
            pairs = [(keep[i + k], v[i]) for i in range(len(v)) if 0 <= i + k < len(keep)]
            agree = sum(a == b for a, b in pairs)
            if agree > best[1]: best = (k, agree, len(pairs))
        rows.append(dict(cnt=cnt_s, veta_session=t, cnt_codes_after_dedup=len(keep), code_vs_veta_offset0=f"{sum(a == b for a, b in zip(keep, v))}/{min(len(keep), len(v))}",
                         best_offset=best[0], code_vs_veta_best=f"{best[1]}/{best[2]}"))
D = pd.DataFrame(rows); D.to_csv(Q / "label_mapping_check_cnt_sessions.csv", index=False, encoding="utf-8-sig"); print(D.to_string(index=False))
