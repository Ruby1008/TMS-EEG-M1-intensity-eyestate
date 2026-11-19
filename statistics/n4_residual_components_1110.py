"""N4 (user 2026-10-28): 11-25 ms residual amplitude vs the ERM components at the group level (17 participants x 6 conditions = 102 cells).
residual = Preproc/QC/QC_summary_wide.csv residual11to25_<lv>_uV (max |ROI average| 11-25 ms; the stage-D metric).
components = ERM ROI P30 27-47, N45 39-59, P60 50-70, N100 120-140 ms (erm_roi, as tep_stats.py) and GMFP 15-60 ms (erm_allch, 61 scalp channels,
SD across channels, as tep_supp_measures.py).
(a) residual ~ eye x intensity rmANOVA (GG) + mean and SD per intensity.
(b) LMMs in R (n4_lmm_1110.R; lme4/lmerTest, ML, Satterthwaite): value ~ eye * intensity (+ residual) + (1 | participant).
(c) repeated-measures correlation residual vs each component (pingouin.rm_corr).
Outputs: _check/N4_cells.csv (input of the R script), _check/N4_anova_residual.csv, _check/N4_rmcorr.csv; N4/N4_residual_vs_components.csv is
written by n4_lmm_1110.R (merges these)."""
import numpy as np, pandas as pd, pingouin as pg
import r1110_lib as L
io = L.io
Q = io.Q[(io.Q.stim == "Real") & (io.Q.status == "processed")].set_index("id")
rows = []
for s in L.S17:
    for e in L.EYES:
        M = L.erm(s, e); te = np.round(M["time"]).astype(int); A = np.asarray(M["erm_allch"], float)    # chan x 3 x time
        for k, lv in enumerate(L.LEV):
            roi = np.asarray(M["erm_roi"][k], float); r = dict(participant=f"{s:03d}", eye=e, intensity=lv, residual=Q.loc[f"spTMS_{s:03d}_Real_{e}", f"residual11to25_{lv}_uV"])
            for c, (a0, a1) in L.WIN_ERM.items(): r[c] = L.wmean(te, roi, a0, a1)
            w = (te >= 15) & (te <= 60); r["GMFP"] = A[:, k, w].std(0).mean(); rows.append(r)
C = pd.DataFrame(rows); assert len(C) == 102 and A.shape[0] == 61
C.to_csv(L.CK / "N4_cells.csv", index=False)
# check against the main component values (ERM P30/N100 of data/Fig6_TEP_MEP_input.csv are the same windows)
TS = pd.read_csv(L.DATA / "TEP_supp_measures_long.csv")
g = TS[(TS.data == "ERM") & (TS.measure == "gmfp_15_60") & TS["sub"].isin(L.S17)].set_index(["sub", "eye", "intensity"]).value
mine = C.assign(sub=C.participant.astype(int)).set_index(["sub", "eye", "intensity"]).GMFP
assert len(g) == 102 and np.allclose(mine.loc[g.index], g), "GMFP differs from data/TEP_supp_measures_long.csv"
a = L.anova2(C.rename(columns={"participant": "sub"}), dv="residual"); a.to_csv(L.CK / "N4_anova_residual.csv", index=False)
desc = C.groupby("intensity").residual.agg(["mean", "std", "median", "min", "max"]); desc.to_csv(L.CK / "N4_residual_descriptives.csv")
rm = []
for c in ("P30", "N45", "P60", "N100", "GMFP"):
    r = pg.rm_corr(data=C, x="residual", y=c, subject="participant"); rm.append(dict(component=c, r_rm=r["r"].iloc[0], df_rm=r["dof"].iloc[0], p_rm=r["pval"].iloc[0]))
pd.DataFrame(rm).to_csv(L.CK / "N4_rmcorr.csv", index=False)
print(a.round(4).to_string()); print(desc.round(2).to_string()); print(pd.DataFrame(rm).round(4).to_string())
