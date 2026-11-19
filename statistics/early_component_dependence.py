"""Dependence of the early ERM components (user 2026-10-19, section 4-B1): cross-participant Spearman correlation matrices between the fixed-window P30
(27-47 ms), N45 (39-59 ms), P60 (50-70 ms) and the 15-60 ms GMFP (data/TEP_supp_measures_long.csv, data = ERM, 17 main participants), per intensity with the
two eye states averaged per participant (main table) and per eye state x intensity (full table). Uncorrected p; n = 17.
Output: data/Early_component_spearman.csv"""
import numpy as np, pandas as pd
from scipy import stats
import tep_io as io

M = pd.read_csv(io.DATA / "TEP_supp_measures_long.csv"); M = M[(M.data == "ERM") & (M.main17) & M.measure.isin(["fixed_P30", "fixed_N45", "fixed_P60", "gmfp_15_60"])]
NAMES = {"fixed_P30": "P30", "fixed_N45": "N45", "fixed_P60": "P60", "gmfp_15_60": "GMFP 15-60"}; order = list(NAMES)
rows = []
def mat(W, lv, eye):
    for i, a in enumerate(order):
        for b in order[i + 1:]:
            r = stats.spearmanr(W[a], W[b]); rows.append(dict(intensity=lv, eye=eye, n=len(W), measure_a=NAMES[a], measure_b=NAMES[b], rho=r.statistic, p=r.pvalue))
for lv in io.LEV:
    g = M[M.intensity == lv]
    W = g.groupby(["sub", "measure"]).value.mean().unstack(); assert len(W) == 17; mat(W, lv, "EO+EC mean")
    for e in io.EYES:
        W = g[g.eye == e].pivot(index="sub", columns="measure", values="value"); assert len(W) == 17; mat(W, lv, e)
R = pd.DataFrame(rows); R.to_csv(io.DATA / "Early_component_spearman.csv", index=False)
pd.set_option("display.width", 200); print(R[R.eye == "EO+EC mean"].round(3).to_string(index=False))
print("per-eye rows:", len(R[R.eye != "EO+EC mean"]), "| range rho (eye-averaged):", R[R.eye == "EO+EC mean"].rho.min().round(2), R[R.eye == "EO+EC mean"].rho.max().round(2))
