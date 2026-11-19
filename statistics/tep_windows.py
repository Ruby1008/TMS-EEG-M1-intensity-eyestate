"""Component windows (decisions 2026-10-16): grand average of the 17 main participants x 6 real conditions, ROI FC3/FC1/C3/C1.
P30: maximum in 20-40 ms (Beck et al., 2024, PMID 39460649); N45: local minimum in 38-55 ms (Belardinelli et al., 2021, PMID 33854132; under the
average reference a relative negativity between P30 and P60, positive amplitude); P60: maximum in 56-82 ms (Belardinelli 2021); window =
peak +-10 ms; N100 fixed 120-140 ms. ERM and _F_TEP (main); BP = _F_ERM band-passed 4.2-36.4 Hz (S-E control, same rule).
Edge check: a peak on a range boundary or not a local extremum is flagged (stop rule for ERM/_F_TEP).
Output: data/TEP_component_windows.csv
"""
import numpy as np, pandas as pd
import tep_io as io

Z = np.load(io.DATA / "ALLCH_avg.npz"); labs = list(Z["labels"]); t = np.round(Z["time"]).astype(int); ri = [labs.index(c) for c in io.ROI]
def ga(kind):
    if kind == "ERM":
        R = io.load_erm(io.MAIN17[0], "EO"); te = np.round(R["time"]).astype(int)
        return te, np.mean([io.load_erm(s, e)["erm_roi"].mean(0) for s in io.MAIN17 for e in io.EYES], 0)
    return t, np.mean([Z[f"spTMS_{s:03d}_Real_{e}|{kind}"][ri].mean(0).mean(1) for s in io.MAIN17 for e in io.EYES], 0)
RANGES = {"P30": (20, 40, 1, "Beck et al., 2024 (PMID 39460649)"), "N45": (38, 55, -1, "Belardinelli et al., 2021 (PMID 33854132)"),
          "P60": (56, 82, 1, "Belardinelli et al., 2021 (PMID 33854132)")}
rows = []
for kind in ("ERM", "TEP", "BP"):
    x, y = ga(kind)
    for c, (a, b, sg, ref) in RANGES.items():
        w = np.where((x >= a) & (x <= b))[0]; i = w[np.argmax(sg * y[w])]; pk = int(x[i])
        loc = sg * y[i] >= sg * y[i - 1] and sg * y[i] >= sg * y[i + 1]
        rows.append(dict(data=kind, component=c, range_ms=f"{a}-{b}", source=ref, peak_ms=pk, peak_uV=float(y[i]), win_start_ms=pk - 10, win_end_ms=pk + 10,
                         edge=pk in (a, b), local_extremum=bool(loc)))
    rows.append(dict(data=kind, component="N100", range_ms="fixed", source="original submission (TW13)", peak_ms=np.nan, peak_uV=float(y[(x >= 120) & (x <= 140)].mean()),
                     win_start_ms=120, win_end_ms=140, edge=False, local_extremum=True))
W = pd.DataFrame(rows); W.to_csv(io.DATA / "TEP_component_windows.csv", index=False); print(W.round(3).to_string())
