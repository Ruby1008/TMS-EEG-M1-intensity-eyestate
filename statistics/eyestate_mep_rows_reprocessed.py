"""Part 3 step 2 (user 2026-10-22): log MEP rows of Supplementary Table S25 (eye-state analyses) and S28 (individual differences) with the reprocessed MEP.
Same code as MF/code/eyestate_null_round8.py (sections 3c (i), (iv), 3f) and eyestate_null_round8.R (3c (ii), (v)), restricted to the measure 'log MEP';
only the MEP input and the output names differ. The N100 and IMF6 ITPC rows of the individual-difference analysis are taken unchanged from
data/EyeState_individual_differences_v2.csv, and Holm over the three tests is recomputed.
Check: with MEPSRC=old the script uses data/MEP_matrix_new.csv and must reproduce the 10-27 table values (log MEP: (i) 2.36/.140, 2.06/.160; (ii) 4.51/.036,
0.65/.524; (iv) 2.76/.131, 1.77/.214; (v) 5.82/.018, 1.38/.243; 3f rho -.04, p .866).
Output: <MEPSET>/analysis/EyeState_MEP_rows_reprocessed.csv, EyeState_individual_differences_reprocessed.csv (or *_oldcheck.csv with MEPSRC=old)"""
data_dir = r"path/to/data"   # release: set this path
import os, subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd, pingouin as pg
from scipy import stats
MF = Path(data_dir + r"\Revision_2026\Manuscript_Figures"); sys.path.insert(0, str(MF / "code"))
import tep_io as io
D = io.DATA; S17 = io.MAIN17; SET = os.environ.get("MEPSET", "main"); OLD = os.environ.get("MEPSRC", "") == "old"
OUT = MF / "MEP_reprocess_2026-10-22" / SET / "analysis"; tag = "_oldcheck" if OLD else "_reprocessed"
MEPm = (pd.read_csv(D / "MEP_matrix_new.csv") if OLD else pd.read_csv(MF / "MEP_reprocess_2026-10-22" / SET / "MEP_matrix_reprocessed.csv", dtype={"subject": str}).assign(sub=lambda x: x.subject.astype(int)).drop(columns="subject")).set_index("sub")
long = []
for s in MEPm.index:
    for e in io.EYES:
        for lv in io.LEV:
            v = MEPm.loc[s, f"{e}_{lv}"]
            if pd.notna(v): long.append(dict(sub=s, eye=e, intensity=lv, measure="log MEP", value=np.log(v)))
L = pd.DataFrame(long); L["main17"] = L["sub"].isin(S17); L.to_csv(OUT / f"_eyestate_mep_long{tag}.csv", index=False)
_q = io.Q[(io.Q.status == "processed") & (io.Q.stim == "Real") & (io.Q.stageD == "sensitivity")]; SENS = [s for s in S17 if s not in set(_q["sub"])]
MEPC = list(MEPm.dropna().index)
ROWS = []
def rm(df, dv, within=("eye", "intensity")):
    return pg.rm_anova(data=df, dv=dv, within=list(within), subject="sub", correction=True, effsize="np2").set_index("Source")
def gg(a, src): r = a.loc[src]; return r["F"], r["ddof1"], r["ddof2"], (r["p_GG_corr"] if "p_GG_corr" in r and pd.notna(r.get("p_GG_corr")) else r["p_unc"]), r["np2"]
for label, subs in (("(i) main analysis", MEPC), ("(iv) sensitivity, no 20-50 uV residual", [s for s in MEPC if s in SENS])):
    g = L[L["sub"].isin(subs)]; full = g.groupby("sub").size(); g = g[g["sub"].isin(full[full == 6].index)]; a = rm(g, "value")
    for src in ("eye", "eye * intensity"):
        F, d1, d2, p, np2 = gg(a, src); ROWS.append(dict(analysis=label, measure="log MEP", effect=src, n=g["sub"].nunique(), F=F, df1=d1, df2=d2, p=p, np2=np2))
Rsrc = r'''suppressMessages({library(lme4); library(lmerTest); library(dplyr)})
L <- read.csv("%s") %%>%% mutate(subject = factor(sub), state = factor(eye, levels = c("EO", "EC")), intens = factor(intensity))
ORD <- read.csv("%s"); L <- L %%>%% left_join(ORD[, c("sub", "order")], by = "sub")
contrasts(L$state) <- contr.sum(2); contrasts(L$intens) <- contr.sum(3); ctl <- lmerControl(optimizer = "bobyqa"); out <- data.frame()
d <- L
f <- lmer(value ~ state * intens + (1 | subject), d, REML = FALSE, control = ctl); a <- as.data.frame(anova(f, type = 3))
for (tm in c("state", "state:intens")) out <- rbind(out, data.frame(analysis = "(ii) LMM, all available participants", measure = "log MEP", effect = tm, n = length(unique(d$subject)), F = a[tm, "F value"], df1 = a[tm, "NumDF"], df2 = a[tm, "DenDF"], p = a[tm, "Pr(>F)"]))
dm <- d %%>%% group_by(sub) %%>%% filter(n() == 6) %%>%% ungroup(); dm$order <- factor(dm$order); contrasts(dm$order) <- contr.sum(2)
f2 <- lmer(value ~ state * intens * order + (1 | subject), dm, REML = FALSE, control = ctl); a2 <- as.data.frame(anova(f2, type = 3))
for (tm in c("state", "state:order")) out <- rbind(out, data.frame(analysis = "(v) block order as between-participant factor (LMM, main set)", measure = "log MEP", effect = tm, n = length(unique(dm$subject)), F = a2[tm, "F value"], df1 = a2[tm, "NumDF"], df2 = a2[tm, "DenDF"], p = a2[tm, "Pr(>F)"]))
write.csv(out, "%s", row.names = FALSE)
''' % ((OUT / f"_eyestate_mep_long{tag}.csv").as_posix(), (D / "EyeState_block_order.csv").as_posix(), (OUT / f"_eyestate_mep_lmm{tag}.csv").as_posix())
(MF / "MEP_reprocess_2026-10-22" / "code" / "_eyestate_mep_rows.R").write_text(Rsrc, encoding="utf-8")
subprocess.run(["Rscript", str(MF / "MEP_reprocess_2026-10-22" / "code" / "_eyestate_mep_rows.R")], check=True)
RES = pd.concat([pd.DataFrame(ROWS), pd.read_csv(OUT / f"_eyestate_mep_lmm{tag}.csv")], ignore_index=True)
RES.to_csv(OUT / f"EyeState_MEP_rows{tag}.csv", index=False)
# 3f: individual differences, log MEP row recomputed; N100 / ITPC rows from the v2 file
LAPd = pd.read_csv(D / "EyeState_laplacian_values.csv")
lap = LAPd[(LAPd.stim == "Real") & (LAPd.version == "after Laplacian (CSD)") & (LAPd.region == "C3")].pivot(index="sub", columns="eye", values="log10_power"); dmu = lap.EC - lap.EO
g = L[L["sub"].isin(S17)].pivot_table(index=["sub", "intensity"], columns="eye", values="value"); y = (g.EC - g.EO).groupby("sub").mean()
j = dmu.index.intersection(y.dropna().index); r = stats.spearmanr(dmu.loc[j], y.loc[j])
V2 = pd.read_csv(D / "EyeState_individual_differences_v2.csv")
ID = pd.concat([pd.DataFrame([dict(x="EC - EO Laplacian mu (C3)", y="EC - EO log MEP", n=len(j), rho=r.statistic, p=r.pvalue)]), V2[~V2.y.str.contains("MEP")][["x", "y", "n", "rho", "p"]]], ignore_index=True)
ID["p_holm"] = pg.multicomp(ID.p.values, method="holm")[1]; ID.to_csv(OUT / f"EyeState_individual_differences{tag}.csv", index=False)
pd.set_option("display.width", 200); print(RES.round(4).to_string(index=False)); print(ID.round(4).to_string(index=False))
