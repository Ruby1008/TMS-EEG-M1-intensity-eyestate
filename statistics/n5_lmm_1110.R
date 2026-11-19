data_dir <- "path/to/data"   # release: set this path
# N5 (user 2026-10-28): single-trial outcome ~ x_within + (1 | subject) per intensity, as section 3e of code/eyestate_null_round8.R (ML, bobyqa,
# all epochs, x centred within participant inside each intensity; change per SD = b x SD of x_within). Input _check/N5_single_trial.csv.
# Variants: (a) P30_main ~ log_mu (must reproduce Table S28: -1.083 / -1.430 / -1.452), (b) P30_IMF45 ~ log_mu, (c) P30_IMF6 ~ log_mu,
# (d) P30_convTEP ~ log_mu, (e) P30_main ~ env_IMF6 (IMF6 pre-stimulus Hilbert envelope, uV) and (e2) P30_main ~ log10 envelope.
suppressMessages({library(lme4); library(lmerTest); library(dplyr)})
RV <- paste0(data_dir, "/Revision_2026/Manuscript_Figures/Review_2026-11-10")
d <- read.csv(file.path(RV, "_check", "N5_single_trial.csv")); d$subject <- factor(d$sub)
ctl <- lmerControl(optimizer = "bobyqa")
V <- list(c("(a) ERM IMF4-7, 27-47 ms (main)", "P30_main", "log_mu"), c("(b) ERM IMF4-5 only, 27-47 ms", "P30_IMF45", "log_mu"),
          c("(c) ERM IMF6 only, 27-47 ms", "P30_IMF6", "log_mu"), c("(d) conventional TEP 0.1-30 Hz, 28-48 ms", "P30_convTEP", "log_mu"),
          c("(e) main P30, predictor IMF6 pre-stimulus envelope (uV)", "P30_main", "env_IMF6"), c("(e2) main P30, predictor log10 IMF6 envelope", "P30_main", "log_env_IMF6"))
out <- data.frame()
for (v in V) for (lv in c(80, 100, 120)) {
  u <- d %>% filter(intensity == lv) %>% group_by(subject) %>% mutate(xw = .data[[v[3]]] - mean(.data[[v[3]]])) %>% ungroup()
  u$y <- u[[v[2]]]
  g <- lmer(y ~ xw + (1 | subject), u, REML = FALSE, control = ctl); cg <- summary(g)$coefficients["xw", ]; sdv <- sd(u$xw)
  out <- rbind(out, data.frame(variant = v[1], outcome = v[2], predictor = v[3], intensity = lv, n = length(unique(u$subject)), n_trials = nrow(u),
                               b = cg["Estimate"], SE = cg["Std. Error"], df = cg["df"], t = cg["t value"], p = cg["Pr(>|t|)"], sd_x_within = sdv, change_per_SD_uV = cg["Estimate"] * sdv,
                               singular = isSingular(g)))
}
write.csv(out, file.path(RV, "_check", "N5_lmm.csv"), row.names = FALSE); print(out, digits = 4)
