data_dir <- "path/to/data"   # release: set this path
# N4 (b) (user 2026-10-28): value ~ eye * intensity (+ residual) + (1 | participant), lme4/lmerTest, ML (REML = FALSE), Satterthwaite df,
# type III F with sum-to-zero contrasts. Input _check/N4_cells.csv (n4_residual_components_1110.py); merges _check/N4_anova_residual.csv and
# _check/N4_rmcorr.csv. Output N4/N4_residual_vs_components.csv.
suppressPackageStartupMessages({library(lme4); library(lmerTest)})
options(contrasts = c("contr.sum", "contr.poly"))
RV <- paste0(data_dir, "/Revision_2026/Manuscript_Figures/Review_2026-11-10")
d <- read.csv(file.path(RV, "_check", "N4_cells.csv"), colClasses = c(participant = "character"))
d$eye <- factor(d$eye, levels = c("EO", "EC")); d$intensity <- factor(d$intensity, levels = c(80, 100, 120))
rm_ <- read.csv(file.path(RV, "_check", "N4_rmcorr.csv"))
out <- list()
for (comp in c("P30", "N45", "P60", "N100", "GMFP")) {
  d$value <- d[[comp]]
  m0 <- lmer(value ~ eye * intensity + (1 | participant), data = d, REML = FALSE)
  m1 <- lmer(value ~ eye * intensity + residual + (1 | participant), data = d, REML = FALSE)
  lrt <- anova(m0, m1)
  rr <- rm_[rm_$component == comp, ]
  for (mn in c("without residual", "with residual")) {
    m <- if (mn == "with residual") m1 else m0
    a <- anova(m, ddf = "Satterthwaite", type = 3)
    co <- summary(m)$coefficients
    for (eff in rownames(a)) {
      out[[length(out) + 1]] <- data.frame(component = comp, model = mn, effect = eff, F = a[eff, "F value"], df1 = a[eff, "NumDF"], df2 = a[eff, "DenDF"],
        p = a[eff, "Pr(>F)"], b_residual = if (mn == "with residual") co["residual", "Estimate"] else NA,
        SE_residual = if (mn == "with residual") co["residual", "Std. Error"] else NA,
        p_residual = if (mn == "with residual") co["residual", "Pr(>|t|)"] else NA,
        LRT_chisq = if (mn == "with residual") lrt$Chisq[2] else NA, LRT_p = if (mn == "with residual") lrt$`Pr(>Chisq)`[2] else NA,
        r_rm = rr$r_rm, df_rm = rr$df_rm, p_rm = rr$p_rm, singular = isSingular(m))
    }
  }
}
res <- do.call(rbind, out)
an <- read.csv(file.path(RV, "_check", "N4_anova_residual.csv"))
res <- rbind(res, data.frame(component = "residual (11-25 ms)", model = "2 x 3 rmANOVA (GG)", effect = an$effect, F = an$F, df1 = an$df1, df2 = an$df2, p = an$p,
  b_residual = NA, SE_residual = NA, p_residual = NA, LRT_chisq = NA, LRT_p = NA, r_rm = NA, df_rm = NA, p_rm = NA, singular = NA))
write.csv(res, file.path(RV, "N4", "N4_residual_vs_components.csv"), row.names = FALSE)
print(res[res$effect %in% c("intensity", "residual", "eye"), c("component", "model", "effect", "F", "df1", "df2", "p", "b_residual", "p_residual", "LRT_p", "r_rm", "p_rm", "singular")], digits = 4)
