data_dir = r"path/to/data"   # release: set this path
# copy made by MEP_reprocess_2026-10-22/code/make_reprocessed_copies.py (2026-10-22); input = reprocessed MEP data (MEPSET = main | excl_artifact); outputs in MEP_reprocess_2026-10-22/<MEPSET>/analysis/.
import sys; sys.path.insert(0, data_dir + r'\Revision_2026\Manuscript_Figures\code')
"""MEP statistics, version A (original submission) versus version B (log) side by side (user 2026-10-19, revised decision 1-f). Not written into the docx.
Data: data/MEP_matrix_new.csv, complete cases (22 participants x 6 cells; participant means of the accepted trials, zero-amplitude trials included).
A: Shapiro-Wilk per cell (raw mV); 2 x 3 rmANOVA (Greenhouse-Geisser) on raw mV; Friedman over the three intensities (eye states averaged); Wilcoxon signed-rank
   post hoc (eye states averaged), Bonferroni over 3 contrasts, rank-biserial r (as mep_posthoc.py).
B: check for zero cell means (listed; no constant is added); 2 x 3 rmANOVA (GG) on ln(MEP); paired t post hoc (eye states averaged), Holm, Cohen's dz;
   Shapiro-Wilk per cell on raw and on log values; residual Shapiro-Wilk (additive participant + cell model) for both scales.
Output: data/MEP_stats_A_vs_B.csv (long table: version, block, item, statistic columns)"""
import numpy as np, pandas as pd, pingouin as pg
from scipy import stats
import tep_io as io
import os as _os
from pathlib import Path as _P
_SET = _os.environ.get("MEPSET", "main"); _RD = _P(data_dir + r"\Revision_2026\Manuscript_Figures\MEP_reprocess_2026-10-22") / _SET
_OUT = _RD / "analysis"; _OUT.mkdir(exist_ok=True)

M = pd.read_csv(_RD / "MEP_matrix_reprocessed.csv", dtype={"subject": str}).rename(columns={"subject": "sub"}).dropna(); N = len(M)
cols = ["EO_80", "EO_100", "EO_120", "EC_80", "EC_100", "EC_120"]
zeros = [(int(r["sub"]), c) for _, r in M.iterrows() for c in cols if r[c] == 0]
print("zero cell means:", zeros if zeros else "none (minimum %.4f mV)" % M[cols].min().min())
rows = []
def add(version, block, item, **kw): rows.append(dict(version=version, block=block, item=item, **kw))
add("B", "zero check", "cells with mean MEP = 0", n=len(zeros), note="none; minimum cell mean %.4f mV" % M[cols].min().min() if not zeros else "; ".join(f"{s} {c}" for s, c in zeros))
def anova(X, label, version):
    long = X.melt(id_vars="sub", var_name="cond", value_name="y"); long["eye"] = long.cond.str.split("_").str[0]; long["intensity"] = long.cond.str.split("_").str[1].astype(int)
    a = pg.rm_anova(data=long, dv="y", within=["eye", "intensity"], subject="sub", correction=True, effsize="np2").set_index("Source")
    for src in ("eye", "intensity", "eye * intensity"):
        r = a.loc[src]; add(version, f"rmANOVA ({label})", src, n=N, F=r["F"], df1=r["ddof1"], df2=r["ddof2"], df1_GG=r["ddof1"] * r["eps"], df2_GG=r["ddof2"] * r["eps"], p_GG=r["p_GG_corr"], p_unc=r["p_unc"], np2=r["np2"], eps=r["eps"])
def sw(X, label, version):
    for c in cols:
        r = stats.shapiro(X[c]); add(version, f"Shapiro–Wilk per cell ({label})", c, n=N, W=r.statistic, p=r.pvalue)
    Y = X[cols].values; res = Y - Y.mean(1, keepdims=True) - Y.mean(0, keepdims=True) + Y.mean(); r = stats.shapiro(res.ravel())
    add(version, f"Shapiro–Wilk residuals ({label})", "participant + cell model residuals", n=res.size, W=r.statistic, p=r.pvalue)
# ---- A
sw(M, "raw mV", "A"); anova(M[["sub"] + cols], "raw mV", "A")
c = {lv: (M[f"EO_{lv}"] + M[f"EC_{lv}"]) / 2 for lv in (80, 100, 120)}
f = stats.friedmanchisquare(c[80], c[100], c[120]); add("A", "Friedman (eye states averaged)", "80 vs 100 vs 120", n=N, chi2=f.statistic, df1=2, p_unc=f.pvalue)
for a_, b_ in ((120, 100), (120, 80), (100, 80)):
    d = c[a_] - c[b_]; w = stats.wilcoxon(c[a_], c[b_]); rk = stats.rankdata(np.abs(d[d != 0])); wp = rk[(d[d != 0] > 0).to_numpy()].sum(); r_rb = (wp - (rk.sum() - wp)) / rk.sum()
    add("A", "post hoc Wilcoxon (eye states averaged, Bonferroni x3)", f"{a_} vs {b_}", n=N, W=w.statistic, p_unc=w.pvalue, p_corr=min(1, 3 * w.pvalue), effect=r_rb, effect_name="rank-biserial r", n_increase=int((d > 0).sum()))
# ---- B
L = M[["sub"] + cols].copy(); L[cols] = np.log(L[cols])
sw(M, "raw mV", "B"); sw(L, "ln mV", "B"); anova(L, "ln mV", "B")
cl = {lv: (L[f"EO_{lv}"] + L[f"EC_{lv}"]) / 2 for lv in (80, 100, 120)}; ph = []
for a_, b_ in ((120, 100), (120, 80), (100, 80)):
    d = cl[a_] - cl[b_]; t = stats.ttest_rel(cl[a_], cl[b_]); ph.append(dict(item=f"{a_} vs {b_}", t=t.statistic, p=t.pvalue, dz=d.mean() / d.std(ddof=1), md=d.mean()))
holm = pg.multicomp([x["p"] for x in ph], method="holm")[1]
for x, ph_ in zip(ph, holm): add("B", "post hoc paired t (eye states averaged, Holm x3)", x["item"], n=N, t=x["t"], df1=N - 1, p_unc=x["p"], p_corr=ph_, effect=x["dz"], effect_name="Cohen's dz", mean_diff_ln=x["md"])
R = pd.DataFrame(rows); R.to_csv(_OUT / "MEP_stats_A_vs_B_reprocessed.csv", index=False, encoding="utf-8-sig")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); print(R.round(4).to_string(index=False))
