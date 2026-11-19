data_dir <- "path/to/data"   # release: set this path
# qc/v2/make_code008.py copy (part 8): spTMS_008_Real_EC no longer excluded.
# copy made by MEP_reprocess_2026-10-22/code/make_reprocessed_copies.py (2026-10-22); MEPSET = main | excl_artifact; outputs in MEP_reprocess_2026-10-22/<MEPSET>/analysis/.
# changes: all accepted trials in the log models, present = p2p >= 0.05 mV, models (ii) and (iii) added (user 2026-10-22, step 7-4)
# Robustness of the single-trial P30 -> log MEP relation (user 2026-10-20, task 1). Not written into the docx.
# Trials: 17 participants, trial2023 mapping, 008_Real_EC excluded, 022_Real_EO trials with mismatching intensity labels removed, accepted trials;
# MEP = 0 excluded from the log models (logistic model separately). Separate models for 100% and 120% RMT; random intercept unless stated; ML; Satterthwaite p.
# Centring within participant (per intensity, on the trials of the model): P30_between = participant mean, P30_within = P30 - participant mean; covariates
# (residual 11-25 ms, mu power) centred within participant the same way; trial order centred within recording (block).
#   (5) main: logMEP ~ P30_within + P30_between + (1 | subject)
#   (2) + resid_within; (3) + order_c; (2)+(3); (6) + mu_within          (all with P30_within and P30_between)
#   (4) logistic: MEP present ~ P30_within + (1 | subject) (all accepted trials, centring on those trials)
#   (7) P30 from the 4.2-36.4 Hz Butterworth band-pass of _F_ERM (no EMD), model (5)
#   (8) % change of MEP per 1 SD of P30_within: exp(b * SD) - 1 (Wald 95% CI); marginal / conditional R2 (Nakagawa, computed from the variance components)
#   additional: random slope (1 + P30_within | subject); individual-peak window P30; both intensities: P30_within x intensity (+ P30_between x intensity)
#   (1) per participant: Spearman rho (P30, MEP), OLS slope of log MEP on P30, BLUP slope (fixed + random slope of the random-slope model)
# Robust (defined in advance by the user): (5) P30_within p < .05 and same sign at both intensities; (1) >= 13/17 participants positive at both intensities
# (sign of the per-participant Spearman rho; OLS and BLUP signs listed too); (2), (3), (6) P30_within still p < .05 with the same sign at both intensities.
# Outputs: data/Fig6_single_trial_robustness.csv, data/Fig6_single_trial_per_participant.csv, data/Fig6_single_trial_zero_mep_by_participant.csv,
#          data/Fig6_single_trial_robust_checks.csv
suppressMessages({library(lme4); library(lmerTest); library(dplyr)})
setwd(paste0(data_dir, "/Revision_2026/Manuscript_Figures"))
SET <- Sys.getenv("MEPSET", "main"); RD <- file.path(paste0(data_dir, "/Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22"), SET); OUTD <- file.path(RD, "analysis"); dir.create(OUTD, showWarnings = FALSE)

X <- read.csv(file.path(OUTD, "Single_trial_robustness_input_reprocessed.csv"))
A <- X %>% filter(!is.na(mep_int_B), mep_int_B == intensity, mep_accept_B == 1) %>%
  mutate(present = as.integer(mep_present_B == 1), subject = factor(sub)) %>% group_by(id) %>% mutate(order_c = trial_order - mean(trial_order)) %>% ungroup()
wc <- function(d, v) { d %>% group_by(subject) %>% mutate("{v}_between" := mean(.data[[v]]), "{v}_within" := .data[[v]] - mean(.data[[v]])) %>% ungroup() }
prep <- function(d) { for (v in c("P30", "P30_BP", "P30_indiv", "resid_11_25", "log_mu")) d <- wc(d, v); d }
OUT <- data.frame(); ctl <- lmerControl(optimizer = "bobyqa")
tidy <- function(m, lv, model, keep = NULL, family = "gaussian") {
  cf <- as.data.frame(summary(m)$coefficients); cf$term <- rownames(cf); if (!is.null(keep)) cf <- cf[cf$term %in% keep, ]
  pcol <- if (family == "gaussian") "Pr(>|t|)" else "Pr(>|z|)"
  data.frame(intensity = lv, model = model, term = cf$term, b = cf$Estimate, se = cf$`Std. Error`, p = cf[[pcol]], n_trials = nobs(m), n_subjects = ngrps(m)[["subject"]],
             singular = isSingular(m), stringsAsFactors = FALSE)
}
r2 <- function(m) { vf <- var(as.vector(model.matrix(m) %*% fixef(m))); vr <- sum(sapply(VarCorr(m), function(v) v[1, 1])); ve <- sigma(m)^2
  c(marginal = vf / (vf + vr + ve), conditional = (vf + vr) / (vf + vr + ve)) }
PP <- data.frame(); EFF <- data.frame()
for (lv in c(100, 120)) {
  u <- prep(A %>% filter(intensity == lv) %>% mutate(logMEP = log(mep_amp_B)))
  K <- c("P30_within", "P30_between", "resid_11_25_within", "order_c", "log_mu_within")
  m5 <- lmer(logMEP ~ P30_within + P30_between + (1 | subject), u, REML = FALSE, control = ctl); OUT <- rbind(OUT, tidy(m5, lv, "(5) within/between, random intercept (main)", K))
  m2 <- lmer(logMEP ~ P30_within + P30_between + resid_11_25_within + (1 | subject), u, REML = FALSE, control = ctl); OUT <- rbind(OUT, tidy(m2, lv, "(2) + residual 11-25 ms", K))
  m3 <- lmer(logMEP ~ P30_within + P30_between + order_c + (1 | subject), u, REML = FALSE, control = ctl); OUT <- rbind(OUT, tidy(m3, lv, "(3) + trial order", K))
  m23 <- lmer(logMEP ~ P30_within + P30_between + resid_11_25_within + order_c + (1 | subject), u, REML = FALSE, control = ctl); OUT <- rbind(OUT, tidy(m23, lv, "(2)+(3) residual + trial order", K))
  m6 <- lmer(logMEP ~ P30_within + P30_between + log_mu_within + (1 | subject), u, REML = FALSE, control = ctl); OUT <- rbind(OUT, tidy(m6, lv, "(6) + mu power", K))
  m7 <- lmer(logMEP ~ P30_BP_within + P30_BP_between + (1 | subject), u, REML = FALSE, control = ctl); OUT <- rbind(OUT, tidy(m7, lv, "(7) band-pass 4.2-36.4 Hz P30", c("P30_BP_within", "P30_BP_between")))
  mi <- lmer(logMEP ~ P30_indiv_within + P30_indiv_between + (1 | subject), u, REML = FALSE, control = ctl); OUT <- rbind(OUT, tidy(mi, lv, "additional: individual-peak window P30", c("P30_indiv_within", "P30_indiv_between")))
  ms <- lmer(logMEP ~ P30_within + P30_between + (1 + P30_within | subject), u, REML = FALSE, control = ctl); OUT <- rbind(OUT, tidy(ms, lv, "additional: random slope", K))
  # (4) logistic
  a <- prep(A %>% filter(intensity == lv))
  g <- glmer(present ~ P30_within + (1 | subject), a, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  OUT <- rbind(OUT, tidy(g, lv, "(4) logistic, MEP present", "P30_within", "binomial"))
  g2 <- glmer(present ~ P30_within + P30_between + (1 | subject), a, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  OUT <- rbind(OUT, tidy(g2, lv, "(iii) logistic, MEP present ~ P30_within + P30_between", c("P30_within", "P30_between"), "binomial"))
  up <- prep(A %>% filter(intensity == lv, present == 1) %>% mutate(logMEP = log(mep_amp_B)))
  mpo <- lmer(logMEP ~ P30_within + P30_between + (1 | subject), up, REML = FALSE, control = ctl); OUT <- rbind(OUT, tidy(mpo, lv, "(ii) trials with an MEP only", K))
  # (8) effect size and R2 for (5)
  b <- fixef(m5)["P30_within"]; se <- sqrt(vcov(m5)["P30_within", "P30_within"]); sdw <- sd(u$P30_within); rr <- r2(m5)
  EFF <- rbind(EFF, data.frame(intensity = lv, sd_P30_within_uV = sdw, pct_change_per_SD = 100 * (exp(b * sdw) - 1), ci_low = 100 * (exp((b - 1.96 * se) * sdw) - 1),
                               ci_high = 100 * (exp((b + 1.96 * se) * sdw) - 1), R2_marginal = rr["marginal"], R2_conditional = rr["conditional"]))
  # (1) per participant
  blup <- coef(ms)$subject
  for (s in levels(droplevels(u$subject))) {
    v <- u %>% filter(subject == s); sp <- suppressWarnings(cor.test(v$P30, v$mep_amp_B, method = "spearman", exact = FALSE)); ol <- summary(lm(logMEP ~ P30, v))$coefficients
    PP <- rbind(PP, data.frame(intensity = lv, sub = as.integer(s), n_trials = nrow(v), rho = unname(sp$estimate), p_rho = sp$p.value, ols_slope = ol["P30", "Estimate"],
                               ols_se = ol["P30", "Std. Error"], ols_p = ol["P30", "Pr(>|t|)"], blup_slope = blup[s, "P30_within"]))
  }
}
# both intensities: P30_within x intensity
u2 <- bind_rows(lapply(c(100, 120), function(lv) prep(A %>% filter(intensity == lv) %>% mutate(logMEP = log(mep_amp_B))))) %>% mutate(Intensity = factor(intensity))
contrasts(u2$Intensity) <- contr.sum(2)
mc <- lmer(logMEP ~ P30_within * Intensity + P30_between * Intensity + (1 | subject), u2, REML = FALSE, control = ctl)
an <- as.data.frame(anova(mc, type = 3)); an$term <- rownames(an)
OUT <- rbind(OUT, data.frame(intensity = "100+120", model = "additional: both intensities, P30_within x intensity (Type III F)", term = an$term, b = an$`F value`, se = NA, p = an$`Pr(>F)`,
                             n_trials = nobs(mc), n_subjects = ngrps(mc)[["subject"]], singular = isSingular(mc)))
write.csv(OUT, file.path(OUTD, "Fig6_single_trial_robustness_reprocessed.csv"), row.names = FALSE); write.csv(PP, file.path(OUTD, "Fig6_single_trial_per_participant_reprocessed.csv"), row.names = FALSE)
write.csv(EFF, file.path(OUTD, "Fig6_single_trial_effect_size_reprocessed.csv"), row.names = FALSE)
Z <- A %>% group_by(intensity, sub, eye) %>% summarise(accepted = n(), mep_zero = sum(present == 0), .groups = "drop"); write.csv(Z, file.path(OUTD, "Fig6_single_trial_no_mep_by_participant_reprocessed.csv"), row.names = FALSE)
# robust checks
w <- function(model) OUT[OUT$model == model & OUT$term %in% c("P30_within"), ]
chk <- function(model) { r <- w(model); all(r$p < .05) && length(unique(sign(r$b))) == 1 }
pos <- PP %>% group_by(intensity) %>% summarise(n = n(), positive_rho = sum(rho > 0), positive_ols = sum(ols_slope > 0), positive_blup = sum(blup_slope > 0), sig_rho = sum(p_rho < .05 & rho > 0),
                                                 sig_neg = sum(p_rho < .05 & rho < 0), median_ols = median(ols_slope), q1_ols = quantile(ols_slope, .25), q3_ols = quantile(ols_slope, .75), .groups = "drop")
C <- data.frame(check = c("(5) P30_within p < .05, same sign, both intensities", "(1) >= 13/17 participants positive (Spearman), both intensities",
                          "(2) + residual 11-25 ms", "(3) + trial order", "(6) + mu power"),
                passed = c(chk("(5) within/between, random intercept (main)"), all(pos$positive_rho >= 13), chk("(2) + residual 11-25 ms"), chk("(3) + trial order"), chk("(6) + mu power")),
                detail = c(paste(sprintf("%s%%: b = %.4f, p = %.4f", w("(5) within/between, random intercept (main)")$intensity, w("(5) within/between, random intercept (main)")$b, w("(5) within/between, random intercept (main)")$p), collapse = "; "),
                           paste(sprintf("%s%%: %d/%d positive (OLS %d, BLUP %d)", pos$intensity, pos$positive_rho, pos$n, pos$positive_ols, pos$positive_blup), collapse = "; "),
                           paste(sprintf("%s%%: b = %.4f, p = %.4f", w("(2) + residual 11-25 ms")$intensity, w("(2) + residual 11-25 ms")$b, w("(2) + residual 11-25 ms")$p), collapse = "; "),
                           paste(sprintf("%s%%: b = %.4f, p = %.4f", w("(3) + trial order")$intensity, w("(3) + trial order")$b, w("(3) + trial order")$p), collapse = "; "),
                           paste(sprintf("%s%%: b = %.4f, p = %.4f", w("(6) + mu power")$intensity, w("(6) + mu power")$b, w("(6) + mu power")$p), collapse = "; ")))
C <- rbind(C, data.frame(check = "robust (all checks passed)", passed = all(C$passed), detail = ""))
write.csv(C, file.path(OUTD, "Fig6_single_trial_robust_checks_reprocessed.csv"), row.names = FALSE); write.csv(pos, file.path(OUTD, "Fig6_single_trial_per_participant_summary_reprocessed.csv"), row.names = FALSE)
print(OUT, digits = 3); print(EFF, digits = 3); print(as.data.frame(pos), digits = 3); print(C)
