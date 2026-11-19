"""AEP-ICA robustness (decision 2026-10-16): P30/N45/P60/N100 before vs after AEP removal, 17 participants, component windows of the main
analysis (data/TEP_component_windows.csv; ERM windows for the ERM, _F_TEP windows for the conventional TEP); 2 x 3 rmANOVA (GG) before and after.
ERM: HHSA (before) / HHSA_AEPr (after); conventional TEP: data/AEP_TEP_roi_waveforms.csv (aep_compare_tep.py).
Output: data/AEP_component_anova.csv, data/AEP_component_values.csv"""
import numpy as np, pandas as pd, pingouin as pg
import tep_io as io
CW = pd.read_csv(io.DATA / "TEP_component_windows.csv"); T = pd.read_csv(io.DATA / "AEP_TEP_roi_waveforms.csv")
def anova(df):
    a = pg.rm_anova(data=df, dv="y", within=["eye", "intensity"], subject="sub", correction=True, effsize="np2").set_index("Source"); out = {}
    for src, k in (("eye", "eye"), ("intensity", "intensity"), ("eye * intensity", "interaction")):
        r = a.loc[src]; out.update({f"F_{k}": r["F"], f"p_{k}": r["p_GG_corr"] if pd.notna(r.get("p_GG_corr", np.nan)) else r["p_unc"], f"np2_{k}": r["np2"], f"eps_{k}": r["eps"] if pd.notna(r.get("eps", np.nan)) else 1.0})
    return out
A, V = [], []
for kind in ("ERM", "TEP"):
    for r in CW[CW.data == kind].itertuples():
        a0, a1 = r.win_start_ms, r.win_end_ms
        for ver in ("before", "after"):
            rows = []
            for s in io.MAIN17:
                for e in io.EYES:
                    if kind == "ERM":
                        R = io.load_erm(s, e, io.HHSA if ver == "before" else io.ROOT / "HHSA_AEPr"); t = np.round(R["time"])
                        for k, lv in enumerate(io.LEV): rows.append(dict(sub=s, eye=e, intensity=lv, y=R["erm_roi"][k][(t >= a0) & (t <= a1)].mean()))
                    else:
                        g = T[(T["sub"] == s) & (T.eye == e) & (T.version == ver)]
                        for lv in io.LEV:
                            h = g[g.intensity == lv]; rows.append(dict(sub=s, eye=e, intensity=lv, y=h[(h.time_ms >= a0) & (h.time_ms <= a1)].uV.mean()))
            df = pd.DataFrame(rows); V.append(df.assign(data=kind, component=r.component, version=ver))
            A.append(dict(data=kind, component=r.component, window_ms=f"{a0}-{a1}", version=ver, n=len(io.MAIN17), **anova(df)))
A = pd.DataFrame(A); A.to_csv(io.DATA / "AEP_component_anova.csv", index=False); pd.concat(V).to_csv(io.DATA / "AEP_component_values.csv", index=False)
pd.set_option("display.width", 200); print(A[["data", "component", "window_ms", "version", "F_intensity", "p_intensity", "np2_intensity", "p_eye", "p_interaction"]].round(4).to_string())
