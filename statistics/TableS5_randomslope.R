data_dir <- "path/to/data"   # release: set this path
# TableS5_randomslope.R (user 2026-11-09, Part A): by-participant random slopes for the single-trial mixed models of Table S5 and Table S6. New file; no existing file is changed.
# Data, trial selection, centring, ML (REML = FALSE) and lmerControl(optimizer = "bobyqa") copied unchanged from
#   TableS5_fullmodel.R            (Table S5: logMEP ~ P30_within + P30_between [+ resid_11_25_within + order_c + log_mu_within] + (1 | subject); MEPSET = main_qc2_008)
#   mu_3e_reprocessed.R            (Table S6: P30 ~ xw + (1 | subject) at 80/100/120% RMT; logMEP ~ xw + (1 | subject) at 100/120% RMT; xw = log_mu centred within participant)
# Step A2: random-intercept models reproduced and checked against TableS5_fullmodel_results.csv / Mu_single_trial_3e_reprocessed.csv (b, SE, p within 1e-6)
#          and the published values (b, SE to 4 decimals, p to 3 decimals; N = 17; 1328 / 1314 trials). Any mismatch -> stop.
# Step A3: same models with (1 + X | subject). If singular or not converged: refit once with (1 + X || subject) and once with control = lmerControl(optimizer = "bobyqa").
#          The original scripts already use bobyqa, so the bobyqa refit has the same specification as the first random-slope fit (reported as such).
# Step A4: fixed effect of X (b, SE, Satterthwaite df, t, p), random-effect SDs and correlation, LRT random intercept vs random slope (anova, refit = FALSE), AIC, BIC, warnings.
# Output: TableS5_S6_randomslope_results.csv (this folder).
suppressMessages({library(lme4); library(lmerTest); library(dplyr)})
MF <- paste0(data_dir, "/Revision_2026/Manuscript_Figures"); setwd(MF)
R8 <- file.path(MF, "MEP_reprocess_2026-10-22/main_qc2_008/analysis"); CD <- file.path(MF, "MEP_reprocess_2026-10-22/qc/v2/code008")
OUTF <- file.path(CD, "TableS5_S6_randomslope_results.csv")
ctl <- lmerControl(optimizer = "bobyqa")

# ---- data, Table S5 (from TableS5_fullmodel.R)
X <- read.csv(file.path(R8, "Single_trial_robustness_input_reprocessed.csv"))
A <- X %>% filter(!is.na(mep_int_B), mep_int_B == intensity, mep_accept_B == 1) %>%
  mutate(present = as.integer(mep_present_B == 1), subject = factor(sub)) %>% group_by(id) %>% mutate(order_c = trial_order - mean(trial_order)) %>% ungroup()
wc <- function(d, v) { d %>% group_by(subject) %>% mutate("{v}_between" := mean(.data[[v]]), "{v}_within" := .data[[v]] - mean(.data[[v]])) %>% ungroup() }
prep <- function(d) { for (v in c("P30", "resid_11_25", "log_mu")) d <- wc(d, v); d }
# ---- data, Table S6 (from mu_3e_reprocessed.R, MEPSET = main_qc2_008)
pw <- read.csv("data/Prestim_power_single_trial.csv"); st <- read.csv("data/ERM_single_trial.csv"); mp <- read.csv(file.path(R8, "EEG_MEP_trial_map_reprocessed.csv"))
dm <- merge(merge(pw, st[, c("id", "epoch", "P30")], by = c("id", "epoch")), mp[, c("id", "epoch", "sub", "eye", "mep_int_B", "mep_amp_B", "mep_accept_B")], by = c("id", "epoch"))
dm$subject <- factor(dm$sub)
mu_mep <- function(lv) dm %>% filter(intensity == lv, !is.na(mep_int_B), mep_int_B == intensity, mep_accept_B == 1, mep_amp_B > 0) %>% mutate(logMEP = log(mep_amp_B)) %>%
  group_by(subject) %>% mutate(xw = log_mu - mean(log_mu)) %>% ungroup()
mu_p30 <- function(lv) dm %>% filter(intensity == lv) %>% group_by(subject) %>% mutate(xw = log_mu - mean(log_mu)) %>% ungroup()

# ---- model list: fixed part, predictor with the random slope, data
S5F <- c(main = "P30_within + P30_between", full = "P30_within + P30_between + resid_11_25_within + order_c + log_mu_within")
ML <- list()
for (lv in c(100, 120)) {
  u <- prep(A %>% filter(intensity == lv) %>% mutate(logMEP = log(mep_amp_B)))
  ML[[length(ML) + 1]] <- list(table = "S5", model = "(a) logMEP ~ P30_within + P30_between", intensity = lv, y = "logMEP", fixed = S5F[["main"]], x = "P30_within", d = u)
  ML[[length(ML) + 1]] <- list(table = "S5", model = "(b) logMEP ~ P30_within + P30_between + resid_11_25_within + order_c + log_mu_within", intensity = lv, y = "logMEP",
                               fixed = S5F[["full"]], x = "P30_within", d = u)
}
for (lv in c(80, 100, 120)) ML[[length(ML) + 1]] <- list(table = "S6", model = "(c) P30 ~ mu_within", intensity = lv, y = "P30", fixed = "xw", x = "xw", d = mu_p30(lv))
for (lv in c(100, 120)) ML[[length(ML) + 1]] <- list(table = "S6", model = "(c) logMEP ~ mu_within", intensity = lv, y = "logMEP", fixed = "xw", x = "xw", d = mu_mep(lv))

fit <- function(f, d, control) {   # fit and collect all warnings, messages and convergence messages (as TableS5_fullmodel.R, plus messages)
  W <- character(0)
  m <- withCallingHandlers(lmer(f, d, REML = FALSE, control = control),
                           warning = function(w) { W <<- c(W, conditionMessage(w)); invokeRestart("muffleWarning") },
                           message = function(w) { W <<- c(W, trimws(conditionMessage(w))); invokeRestart("muffleMessage") })
  W <- c(W, m@optinfo$conv$lme4$messages, m@optinfo$warnings)
  sing <- isSingular(m); if (sing && !any(grepl("singular", W))) W <- c(W, "singular fit (isSingular = TRUE)")
  conv_ok <- length(m@optinfo$conv$lme4$messages) == 0 && m@optinfo$conv$opt == 0
  list(m = m, warn = if (length(W)) paste(unique(unlist(W)), collapse = " | ") else "none", singular = sing, converged = conv_ok)
}

# ---- step A2: reproduction of the random-intercept models
PUB <- data.frame(model = c("(a)", "(b)", "(a)", "(b)", "mep", "mep"), intensity = c(100, 100, 120, 120, 100, 120),
                  b = c(0.0327, 0.0292, 0.0115, 0.0116, -0.101, -0.028), se = c(0.0069, 0.0071, 0.0045, 0.0045, NA, NA), p = c(NA, NA, 0.011, 0.011, 0.165, 0.675))
F5 <- read.csv(file.path(CD, "TableS5_fullmodel_results.csv")); F6 <- read.csv(file.path(R8, "Mu_single_trial_3e_reprocessed.csv"))
MISM <- character(0); RI <- list()
for (k in seq_along(ML)) {
  s <- ML[[k]]; r <- fit(as.formula(paste(s$y, "~", s$fixed, "+ (1 | subject)")), s$d, ctl); RI[[k]] <- r
  cf <- summary(r$m)$coefficients[s$x, ]; tag <- sprintf("%s %s %d%%", s$table, s$model, s$intensity)
  if (s$table == "S5") {
    o <- F5[F5$intensity == s$intensity & F5$term == "P30_within" & F5$section == ifelse(grepl("^\\(a\\)", s$model), "step2_reproduction", "step3_full_model") &
              F5$model %in% c("(5) within/between, random intercept (main)", "full"), ]
    pb <- PUB[PUB$model == substr(s$model, 1, 3) & PUB$intensity == s$intensity, ]
  } else {
    o <- F6[F6$measure == "mu (FC3/FC1/C3/C1)" & F6$effect == paste0(s$intensity, "% RMT") &
              F6$analysis == ifelse(s$y == "P30", "single-trial P30 ~ power_within + (1 | subject)", "single-trial log MEP ~ power_within + (1 | subject)"), ]
    names(o)[names(o) == "se"] <- "SE"; pb <- if (s$y == "logMEP") PUB[PUB$model == "mep" & PUB$intensity == s$intensity, ] else PUB[0, ]
  }
  if (nrow(o) != 1 || abs(o$b - cf[["Estimate"]]) > 1e-6 || abs(o$SE - cf[["Std. Error"]]) > 1e-6 || abs(o$p - cf[["Pr(>|t|)"]]) > 1e-6)
    MISM <- c(MISM, sprintf("%s: differs from the result file", tag))
  if (nrow(pb) == 1) {
    nd <- ifelse(s$table == "S5", "%.4f", "%.3f")
    if (sprintf(nd, cf[["Estimate"]]) != sprintf(nd, pb$b)) MISM <- c(MISM, sprintf("%s: b %.6f vs published %s", tag, cf[["Estimate"]], pb$b))
    if (!is.na(pb$se) && sprintf("%.4f", cf[["Std. Error"]]) != sprintf("%.4f", pb$se)) MISM <- c(MISM, sprintf("%s: SE %.6f vs published %.4f", tag, cf[["Std. Error"]], pb$se))
    if (is.na(pb$p) && cf[["Pr(>|t|)"]] >= .001) MISM <- c(MISM, sprintf("%s: p %.5f vs published < .001", tag, cf[["Pr(>|t|)"]]))
    if (!is.na(pb$p) && sprintf("%.3f", cf[["Pr(>|t|)"]]) != sprintf("%.3f", pb$p)) MISM <- c(MISM, sprintf("%s: p %.5f vs published %.3f", tag, cf[["Pr(>|t|)"]], pb$p))
  }
  ns <- length(unique(s$d$subject)); if (ns != 17) MISM <- c(MISM, sprintf("%s: %d participants (expected 17)", tag, ns))
  if (s$table == "S5" && nobs(r$m) != c("100" = 1328, "120" = 1314)[[as.character(s$intensity)]]) MISM <- c(MISM, sprintf("%s: %d trials", tag, nobs(r$m)))
}
if (length(MISM)) { cat("STOP: checkpoint A2 failed\n", paste(MISM, collapse = "\n"), "\n"); stop("A2 mismatch; A3-A4 not run") }
cat("A2: all random-intercept models reproduced (result files within 1e-6; published values to the reported decimals; N = 17; 1328 / 1314 trials)\n")

# ---- steps A3-A4
re_sd <- function(m, x) {
  vc <- as.data.frame(VarCorr(m)); g <- vc[vc$grp != "Residual", ]
  sdi <- g$sdcor[is.na(g$var2) & g$var1 == "(Intercept)"]; sds <- g$sdcor[is.na(g$var2) & g$var1 == x]; cr <- g$sdcor[!is.na(g$var2)]
  c(sd_intercept = if (length(sdi)) sdi else NA, sd_slope = if (length(sds)) sds else NA, cor_int_slope = if (length(cr)) cr else NA, sd_residual = vc$sdcor[vc$grp == "Residual"])
}
row <- function(s, variant, r, r0) {
  cf <- summary(r$m)$coefficients[s$x, ]; rs <- re_sd(r$m, s$x)
  lrt <- if (is.null(r0)) NULL else as.data.frame(anova(r0$m, r$m, refit = FALSE))
  data.frame(table = s$table, model = s$model, intensity = s$intensity, predictor = ifelse(s$x == "xw", "mu_within (log_mu centred within participant)", s$x),
             random_effects = variant, n_participants = length(unique(s$d$subject)), n_trials = nobs(r$m),
             b = cf[["Estimate"]], SE = cf[["Std. Error"]], df = cf[["df"]], t = cf[["t value"]], p = cf[["Pr(>|t|)"]],
             b_4dp = sprintf("%.4f", cf[["Estimate"]]), SE_4dp = sprintf("%.4f", cf[["Std. Error"]]),
             p_3dp = ifelse(cf[["Pr(>|t|)"]] < .001, "< .001", sub("^0", "", sprintf("%.3f", cf[["Pr(>|t|)"]]))),
             sd_intercept = rs[["sd_intercept"]], sd_slope = rs[["sd_slope"]], cor_int_slope = rs[["cor_int_slope"]], sd_residual = rs[["sd_residual"]],
             logLik = as.numeric(logLik(r$m)), AIC = AIC(r$m), BIC = BIC(r$m),
             LRT_vs_random_intercept_chisq = if (is.null(lrt)) NA else lrt$Chisq[2], LRT_df = if (is.null(lrt)) NA else lrt$Df[2], LRT_p = if (is.null(lrt)) NA else lrt$`Pr(>Chisq)`[2],
             AIC_random_intercept = if (is.null(r0)) NA else AIC(r0$m), BIC_random_intercept = if (is.null(r0)) NA else BIC(r0$m),
             singular = r$singular, converged = r$converged, warnings = r$warn, stringsAsFactors = FALSE)
}
OUT <- data.frame()
for (k in seq_along(ML)) {
  s <- ML[[k]]; r0 <- RI[[k]]
  OUT <- rbind(OUT, row(s, "(1 | participant) [published]", r0, NULL))
  r1 <- fit(as.formula(sprintf("%s ~ %s + (1 + %s | subject)", s$y, s$fixed, s$x)), s$d, ctl)
  OUT <- rbind(OUT, row(s, "(1 + X | participant), bobyqa [as original scripts]", r1, r0))
  if (r1$singular || !r1$converged) {
    r2 <- fit(as.formula(sprintf("%s ~ %s + (1 + %s || subject)", s$y, s$fixed, s$x)), s$d, ctl)
    OUT <- rbind(OUT, row(s, "refit 1: (1 + X || participant), uncorrelated, bobyqa", r2, r0))
    r3 <- fit(as.formula(sprintf("%s ~ %s + (1 + %s | subject)", s$y, s$fixed, s$x)), s$d, lmerControl(optimizer = "bobyqa"))
    OUT <- rbind(OUT, row(s, "refit 2: (1 + X | participant), control = lmerControl(optimizer = 'bobyqa') [same specification as the first random-slope fit]", r3, r0))
  }
}
write.csv(OUT, OUTF, row.names = FALSE)

# ---- plain-text summary
for (i in seq_len(nrow(OUT))) with(OUT[i, ], cat(sprintf(
  "\n%s %s | %d%% RMT | %s\n  N = %d, trials = %d | b = %s, SE = %s, df = %.1f, t = %.3f, p %s\n  SD int = %.4f, SD slope = %s, r = %s, SD resid = %.4f | AIC = %.2f, BIC = %.2f%s\n  warnings: %s\n",
  table, model, intensity, random_effects, n_participants, n_trials, b_4dp, SE_4dp, df, t, ifelse(p_3dp == "< .001", "< .001", paste("=", p_3dp)),
  sd_intercept, ifelse(is.na(sd_slope), "NA", sprintf("%.4f", sd_slope)), ifelse(is.na(cor_int_slope), "NA", sprintf("%.3f", cor_int_slope)), sd_residual, AIC, BIC,
  ifelse(is.na(LRT_vs_random_intercept_chisq), "", sprintf(" | LRT chi2(%d) = %.3f, p = %.4g", LRT_df, LRT_vs_random_intercept_chisq, LRT_p)), warnings)))
cat("\nwritten:", OUTF, "\n")
