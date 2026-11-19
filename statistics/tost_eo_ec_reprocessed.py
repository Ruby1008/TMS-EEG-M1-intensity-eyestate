data_dir = r"path/to/data"   # release: set this path
# copy made by MEP_reprocess_2026-10-22/code/make_reprocessed_copies.py (2026-10-22); input = reprocessed MEP data (MEPSET = main | excl_artifact); outputs in MEP_reprocess_2026-10-22/<MEPSET>/analysis/.
import sys; sys.path.insert(0, data_dir + r'\Revision_2026\Manuscript_Figures\code')
# v2 copy (make_itpc_v2_copies.py, 2026-10-22): ITPC windows of rule v2 (data/Fig8_ITPC_N100_input_v2.csv); outputs TOST_EO_EC_v2.csv and EyeState_TOST_bounds_v2.csv.
"""Equivalence tests EO vs EC (user 2026-10-18, section 3.2; pre-specified analysis added for the reviewers).
Paired two one-sided tests (TOST; Lakens 2017) on the EO - EC difference, each intensity separately, equivalence bounds dz = ±0.5
(raw bound = 0.5 x SD of the paired difference), 90% CI of the mean difference, TOST p = max of the two one-sided p; Holm over all 21 TOST tests.
Measures: log MEP amplitude (natural log of the participant x condition mean, mV; data/MEP_matrix_new.csv, complete cases n = 22; version B, user 2026-10-20); ERM P30 27-47, N45 39-59, P60 50-70, N100 120-140 ms
(data/TEP_component_values.csv, data = ERM, n = 17); ITPC IMF6 and IMF7, main-analysis windows (longest significant segment; data/Fig8_ITPC_N100_input.csv, n = 17).
The ordinary paired t-test (difference from 0) is listed alongside. Output: data/TOST_EO_EC.csv"""
import numpy as np, pandas as pd
from scipy import stats
import tep_io as io
import os as _os
from pathlib import Path as _P
_SET = _os.environ.get("MEPSET", "main"); _RD = _P(data_dir + r"\Revision_2026\Manuscript_Figures\MEP_reprocess_2026-10-22") / _SET
_OUT = _RD / "analysis"; _OUT.mkdir(exist_ok=True)

D = io.DATA; rows = []
def tost(eo, ec, measure, lv, unit):
    d = np.asarray(eo, float) - np.asarray(ec, float); n = len(d); m = d.mean(); sd = d.std(ddof=1); se = sd / np.sqrt(n); df = n - 1
    lo, hi = -0.5 * sd, 0.5 * sd
    p_lo = stats.t.sf((m - lo) / se, df); p_hi = stats.t.cdf((m - hi) / se, df)
    ci = m + np.array([-1, 1]) * stats.t.ppf(.95, df) * se; t0 = stats.ttest_1samp(d, 0)
    rows.append(dict(measure=measure, intensity=lv, unit=unit, n=n, mean_EO=np.mean(eo), mean_EC=np.mean(ec), mean_diff_EO_minus_EC=m, sd_diff=sd, dz=m / sd,
                     bound_low=lo, bound_high=hi, ci90_low=ci[0], ci90_high=ci[1], t_lower=(m - lo) / se, t_upper=(m - hi) / se, df=df,
                     p_lower=p_lo, p_upper=p_hi, p_tost=max(p_lo, p_hi), t_diff=t0.statistic, p_diff=t0.pvalue))
M = pd.read_csv(_RD / "MEP_matrix_reprocessed.csv", dtype={"subject": str}).rename(columns={"subject": "sub"}).dropna()
for lv in io.LEV: tost(np.log(M[f"EO_{lv}"]), np.log(M[f"EC_{lv}"]), "log MEP", lv, "ln mV")   # 2026-10-20 (user decision 1-4): MEP on the log scale (version B)
CV = pd.read_csv(D / "TEP_component_values.csv"); CV = CV[(CV.data == "ERM") & CV["sub"].isin(io.MAIN17)]
for c in ("P30", "N45", "P60", "N100"):
    for lv in io.LEV:
        w = CV[(CV.component == c) & (CV.intensity == lv)].pivot(index="sub", columns="eye", values="y"); assert len(w) == 17
        tost(w.EO, w.EC, f"ERM {c}", lv, "µV")
F8 = pd.read_csv(D / "Fig8_ITPC_N100_input_v2.csv"); assert len(F8) == 17
for imf in (6, 7):
    for lv in io.LEV: tost(F8[f"ITPC_IMF{imf}_EO_{lv}"], F8[f"ITPC_IMF{imf}_EC_{lv}"], f"ITPC IMF{imf}", lv, "ITPC")
T = pd.DataFrame(rows)
p = T.p_tost.values; o = np.argsort(p); h = np.minimum(1, np.maximum.accumulate(p[o] * (len(p) - np.arange(len(p))))); T.loc[T.index[o], "p_tost_holm"] = h
T["equivalent_holm"] = T.p_tost_holm < .05; T.to_csv(_OUT / "TOST_EO_EC_reprocessed.csv", index=False)
# smallest passing bounds (as eyestate_null_round8.py section 3d)
tb = []
for r in T.itertuples():
    se = r.sd_diff / np.sqrt(r.n); df = r.n - 1
    for alpha, lab in ((.05, ".05"), (.05 / 21, ".05/21")):
        b = abs(r.mean_diff_EO_minus_EC) + stats.t.ppf(1 - alpha, df) * se
        tb.append(dict(measure=r.measure, intensity=r.intensity, unit=r.unit, n=r.n, mean_diff_EO_minus_EC=r.mean_diff_EO_minus_EC, sd_diff=r.sd_diff, dz=r.dz, alpha=lab,
                       original_bound_dz=0.5, original_bound_raw=r.bound_high, p_tost_original=r.p_tost, p_tost_holm_original=r.p_tost_holm, smallest_passing_bound_raw=b, smallest_passing_bound_dz=b / r.sd_diff))
pd.DataFrame(tb).to_csv(_OUT / "EyeState_TOST_bounds_reprocessed.csv", index=False)
print(T[["measure", "intensity", "n", "mean_diff_EO_minus_EC", "dz", "ci90_low", "ci90_high", "bound_high", "p_tost", "p_tost_holm", "p_diff"]].round(4).to_string(index=False))
print("equivalent after Holm:", int(T.equivalent_holm.sum()), "of", len(T), "; uncorrected:", int((T.p_tost < .05).sum()))
