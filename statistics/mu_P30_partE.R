data_dir <- "path/to/data"   # release: set this path
# mu_P30_partE.R (user 2026-11-10, Part E): pre-stimulus mu -> single-trial P30, heterogeneity of the slopes. New file; no existing file is changed.
# E1: per-participant OLS slopes of ERM P30 (IMF4-7, 27-47 ms) on mu_within, trials / centring / exclusions as mu_3e_reprocessed.R (MEPSET = main_qc2_008):
#     all merged epochs of the intensity, log_mu centred within participant inside the intensity. Per participant: b, SE, 95% CI (t, n - 2 df).
#     Per intensity: n negative (of 17), n CI excluding 0 (negative / positive), median and IQR (quantile type 7), one-sample t test and Wilcoxon signed-rank test vs 0.
#     Checkpoint: the random-intercept model on the same data reproduces Mu_single_trial_3e_reprocessed.csv (b, SE, p within 1e-6).
#     Output: mu_P30_per_participant_slopes.csv (this folder).
# E2: outcome-specific P30 models of Review_2026-11-10/code/n5_lmm_1110.R (input Review_2026-11-10/_check/N5_single_trial.csv, all epochs, ML, bobyqa):
#     (b) P30_IMF45 ~ mu_within, (c) P30_IMF6 ~ mu_within (alpha IMF), (d) P30_convTEP ~ mu_within, each + (1 | subject) per intensity.
#     Checkpoint: random-intercept fits reproduce Review_2026-11-10/N5/N5_mu_P30_variants.csv (b, SE, p within 1e-6).
#     Then (1 + mu_within | subject); if singular or not converged: refit once with (1 + mu_within || subject) and once with lmerControl(optimizer = "bobyqa")
#     (the original already uses bobyqa, so that refit has the same specification as the first random-slope fit). LRT via anova(refit = FALSE).
#     Output: mu_P30_IMF_randomslope.csv (this folder).
suppressMessages({library(lme4); library(lmerTest); library(dplyr)})
MF <- paste0(data_dir, "/Revision_2026/Manuscript_Figures"); setwd(MF)
R8 <- file.path(MF, "MEP_reprocess_2026-10-22/main_qc2_008/analysis"); CD <- file.path(MF, "MEP_reprocess_2026-10-22/qc/v2/code008")
RV <- file.path(MF, "Review_2026-11-10")
ctl <- lmerControl(optimizer = "bobyqa")
fmt_p <- function(p) ifelse(p < .001, "< .001", sub("^0", "", sprintf("%.3f", p)))

# ---- E1 data (from mu_3e_reprocessed.R)
pw <- read.csv("data/Prestim_power_single_trial.csv"); st <- read.csv("data/ERM_single_trial.csv"); mp <- read.csv(file.path(R8, "EEG_MEP_trial_map_reprocessed.csv"))
dm <- merge(merge(pw, st[, c("id", "epoch", "P30")], by = c("id", "epoch")), mp[, c("id", "epoch", "sub", "eye", "mep_int_B", "mep_amp_B", "mep_accept_B")], by = c("id", "epoch"))
dm$subject <- factor(dm$sub)
F6 <- read.csv(file.path(R8, "Mu_single_trial_3e_reprocessed.csv"))
P1 <- data.frame(); S1 <- data.frame(); MISM <- character(0)
for (lv in c(80, 100, 120)) {
  v <- dm %>% filter(intensity == lv) %>% group_by(subject) %>% mutate(xw = log_mu - mean(log_mu)) %>% ungroup()
  g <- lmer(P30 ~ xw + (1 | subject), v, REML = FALSE, control = ctl); cg <- summary(g)$coefficients["xw", ]
  o <- F6[F6$measure == "mu (FC3/FC1/C3/C1)" & F6$effect == paste0(lv, "% RMT") & F6$analysis == "single-trial P30 ~ power_within + (1 | subject)", ]
  if (nrow(o) != 1 || abs(o$b - cg[["Estimate"]]) > 1e-6 || abs(o$se - cg[["Std. Error"]]) > 1e-6 || abs(o$p - cg[["Pr(>|t|)"]]) > 1e-6 || nrow(v) != o$n_obs)
    MISM <- c(MISM, sprintf("E1 %d%%: random-intercept model differs from Mu_single_trial_3e_reprocessed.csv", lv))
  for (s in levels(droplevels(v$subject))) {
    w <- v[v$subject == s, ]; f <- lm(P30 ~ xw, w); cf <- summary(f)$coefficients["xw", ]; ci <- confint(f, "xw", level = .95)
    P1 <- rbind(P1, data.frame(section = "E1_participant", intensity = lv, participant = as.integer(s), n_trials = nrow(w), b = cf[["Estimate"]], SE = cf[["Std. Error"]],
                               ci_low = ci[1], ci_high = ci[2], t = cf[["t value"]], df = f$df.residual, p = cf[["Pr(>|t|)"]],
                               ci_excludes_0 = ifelse(ci[2] < 0, "negative", ifelse(ci[1] > 0, "positive", "no"))))
  }
  b <- P1$b[P1$intensity == lv]; ex <- P1$ci_excludes_0[P1$intensity == lv]; tt <- t.test(b, mu = 0); wt <- wilcox.test(b, mu = 0)
  S1 <- rbind(S1, data.frame(section = "E1_summary", intensity = lv, n_participants = length(b), n_trials = nrow(v), n_negative = sum(b < 0),
                             n_ci_excl0_negative = sum(ex == "negative"), n_ci_excl0_positive = sum(ex == "positive"),
                             median = median(b), q25 = unname(quantile(b, .25)), q75 = unname(quantile(b, .75)), IQR = IQR(b), mean = mean(b), sd = sd(b),
                             t = unname(tt$statistic), t_df = unname(tt$parameter), t_p = tt$p.value, wilcoxon_V = unname(wt$statistic), wilcoxon_p = wt$p.value,
                             wilcoxon_exact = length(b) < 50 && !any(duplicated(b)) && !any(b == 0),   # R default: exact p for n < 50 without ties or zeros
                             RI_b = cg[["Estimate"]], RI_SE = cg[["Std. Error"]], RI_p = cg[["Pr(>|t|)"]]))
}

# ---- E2 checkpoint (n5_lmm_1110.R)
d5 <- read.csv(file.path(RV, "_check", "N5_single_trial.csv")); d5$subject <- factor(d5$sub)
N5 <- read.csv(file.path(RV, "N5", "N5_mu_P30_variants.csv"))
VV <- list(c("(b) ERM IMF4-5 only, 27-47 ms", "P30_IMF45"), c("(c) ERM IMF6 only, 27-47 ms", "P30_IMF6"), c("(d) conventional TEP 0.1-30 Hz, 28-48 ms", "P30_convTEP"))
ML <- list()
for (vv in VV) for (lv in c(80, 100, 120)) {
  u <- d5 %>% filter(intensity == lv) %>% group_by(subject) %>% mutate(xw = log_mu - mean(log_mu)) %>% ungroup(); u$y <- u[[vv[2]]]
  ML[[length(ML) + 1]] <- list(variant = vv[1], outcome = vv[2], intensity = lv, d = u)
}
fit <- function(f, d, control) {   # as TableS5_randomslope.R
  W <- character(0)
  m <- withCallingHandlers(lmer(f, d, REML = FALSE, control = control),
                           warning = function(w) { W <<- c(W, conditionMessage(w)); invokeRestart("muffleWarning") },
                           message = function(w) { W <<- c(W, trimws(conditionMessage(w))); invokeRestart("muffleMessage") })
  W <- c(W, m@optinfo$conv$lme4$messages, m@optinfo$warnings)
  sing <- isSingular(m); if (sing && !any(grepl("singular", W))) W <- c(W, "singular fit (isSingular = TRUE)")
  conv_ok <- length(m@optinfo$conv$lme4$messages) == 0 && m@optinfo$conv$opt == 0
  list(m = m, warn = if (length(W)) paste(unique(unlist(W)), collapse = " | ") else "none", singular = sing, converged = conv_ok)
}
RI <- list()
for (k in seq_along(ML)) {
  s <- ML[[k]]; r <- fit(y ~ xw + (1 | subject), s$d, ctl); RI[[k]] <- r; cf <- summary(r$m)$coefficients["xw", ]
  o <- N5[N5$variant == s$variant & N5$intensity == s$intensity, ]
  if (nrow(o) != 1 || abs(o$b - cf[["Estimate"]]) > 1e-6 || abs(o$SE - cf[["Std. Error"]]) > 1e-6 || abs(o$p - cf[["Pr(>|t|)"]]) > 1e-6 || nobs(r$m) != o$n_trials)
    MISM <- c(MISM, sprintf("E2 %s %d%%: differs from N5_mu_P30_variants.csv", s$variant, s$intensity))
}
if (length(MISM)) { cat("STOP: checkpoint failed\n", paste(MISM, collapse = "\n"), "\n"); stop("checkpoint mismatch; nothing written") }
cat("checkpoints passed: E1 random-intercept models = Mu_single_trial_3e_reprocessed.csv; E2 random-intercept models = N5_mu_P30_variants.csv (within 1e-6)\n")
write.csv(bind_rows(P1, S1), file.path(CD, "mu_P30_per_participant_slopes.csv"), row.names = FALSE)

# ---- E2 random slopes
re_sd <- function(m) {
  vc <- as.data.frame(VarCorr(m)); g <- vc[vc$grp != "Residual", ]
  sdi <- g$sdcor[is.na(g$var2) & g$var1 == "(Intercept)"]; sds <- g$sdcor[is.na(g$var2) & g$var1 == "xw"]; cr <- g$sdcor[!is.na(g$var2)]
  c(sd_intercept = if (length(sdi)) sdi else NA, sd_slope = if (length(sds)) sds else NA, cor_int_slope = if (length(cr)) cr else NA, sd_residual = vc$sdcor[vc$grp == "Residual"])
}
row <- function(s, variant, r, r0) {
  cf <- summary(r$m)$coefficients["xw", ]; rs <- re_sd(r$m); lrt <- if (is.null(r0)) NULL else as.data.frame(anova(r0$m, r$m, refit = FALSE))
  data.frame(variant = s$variant, outcome = s$outcome, predictor = "mu_within (log_mu centred within participant)", intensity = s$intensity, random_effects = variant,
             n_participants = length(unique(s$d$subject)), n_trials = nobs(r$m), b = cf[["Estimate"]], SE = cf[["Std. Error"]], df = cf[["df"]], t = cf[["t value"]], p = cf[["Pr(>|t|)"]],
             b_3dp = sprintf("%.3f", cf[["Estimate"]]), SE_3dp = sprintf("%.3f", cf[["Std. Error"]]), p_3dp = fmt_p(cf[["Pr(>|t|)"]]),
             sd_intercept = rs[["sd_intercept"]], sd_slope = rs[["sd_slope"]], cor_int_slope = rs[["cor_int_slope"]], sd_residual = rs[["sd_residual"]],
             logLik = as.numeric(logLik(r$m)), AIC = AIC(r$m), BIC = BIC(r$m),
             LRT_vs_random_intercept_chisq = if (is.null(lrt)) NA else lrt$Chisq[2], LRT_df = if (is.null(lrt)) NA else lrt$Df[2], LRT_p = if (is.null(lrt)) NA else lrt$`Pr(>Chisq)`[2],
             AIC_random_intercept = if (is.null(r0)) NA else AIC(r0$m), BIC_random_intercept = if (is.null(r0)) NA else BIC(r0$m),
             singular = r$singular, converged = r$converged, warnings = r$warn, stringsAsFactors = FALSE)
}
OUT <- data.frame()
for (k in seq_along(ML)) {
  s <- ML[[k]]; r0 <- RI[[k]]
  OUT <- rbind(OUT, row(s, "(1 | participant) [as Review_2026-11-10 N5]", r0, NULL))
  r1 <- fit(y ~ xw + (1 + xw | subject), s$d, ctl); OUT <- rbind(OUT, row(s, "(1 + mu_within | participant), bobyqa [as original]", r1, r0))
  if (r1$singular || !r1$converged) {
    r2 <- fit(y ~ xw + (1 + xw || subject), s$d, ctl); OUT <- rbind(OUT, row(s, "refit 1: (1 + mu_within || participant), uncorrelated, bobyqa", r2, r0))
    r3 <- fit(y ~ xw + (1 + xw | subject), s$d, lmerControl(optimizer = "bobyqa"))
    OUT <- rbind(OUT, row(s, "refit 2: (1 + mu_within | participant), control = lmerControl(optimizer = 'bobyqa') [same specification as the first random-slope fit]", r3, r0))
  }
}
write.csv(OUT, file.path(CD, "mu_P30_IMF_randomslope.csv"), row.names = FALSE)

# ---- plain-text summary
for (i in seq_len(nrow(S1))) with(S1[i, ], cat(sprintf(
  "\nE1 %d%% RMT: N = %d, trials = %d | negative slopes %d/%d | CI excludes 0: %d negative, %d positive | median %.3f, IQR %.3f to %.3f (IQR %.3f) | t(%d) = %.3f, p %s | Wilcoxon V = %d, p %s (exact = %s)\n",
  intensity, n_participants, n_trials, n_negative, n_participants, n_ci_excl0_negative, n_ci_excl0_positive, median, q25, q75, IQR, t_df, t,
  ifelse(t_p < .001, "< .001", paste("=", fmt_p(t_p))), wilcoxon_V, ifelse(wilcoxon_p < .001, "< .001", paste("=", fmt_p(wilcoxon_p))), wilcoxon_exact)))
for (i in seq_len(nrow(OUT))) with(OUT[i, ], cat(sprintf(
  "\nE2 %s | %d%% RMT | %s\n  N = %d, trials = %d | b = %s, SE = %s, df = %.1f, t = %.3f, p %s\n  SD int = %.4f, SD slope = %s, r = %s, SD resid = %.4f | AIC = %.2f, BIC = %.2f%s\n  warnings: %s\n",
  variant, intensity, random_effects, n_participants, n_trials, b_3dp, SE_3dp, df, t, ifelse(p_3dp == "< .001", "< .001", paste("=", p_3dp)),
  sd_intercept, ifelse(is.na(sd_slope), "NA", sprintf("%.4f", sd_slope)), ifelse(is.na(cor_int_slope), "NA", sprintf("%.3f", cor_int_slope)), sd_residual, AIC, BIC,
  ifelse(is.na(LRT_vs_random_intercept_chisq), "", sprintf(" | LRT chi2(%d) = %.3f, p = %.4g", LRT_df, LRT_vs_random_intercept_chisq, LRT_p)), warnings)))
cat("\nwritten:", file.path(CD, "mu_P30_per_participant_slopes.csv"), "and", file.path(CD, "mu_P30_IMF_randomslope.csv"), "\n")
