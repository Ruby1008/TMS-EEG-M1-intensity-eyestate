data_dir <- "path/to/data"   # release: set this path
# TableS5_fullmodel.R (user 2026-10-29): "full model" sensitivity analysis for Table S5 (single-trial P30 -> log MEP). New file; no existing file is changed.
# Data, trial selection, centring and covariates copied from single_trial_robustness_reprocessed.R in this folder (MEPSET = main_qc2_008 = Table S5):
#   accepted trials with matching intensity label; logMEP = ln(mep_amp_B); P30, resid_11_25, log_mu (log10 8-12 Hz) centred within participant per intensity;
#   order_c = trial_order centred within recording (block). ML (REML = FALSE), bobyqa, Satterthwaite df (lmerTest).
# Step 2: models (5) (2) (3) (6) reproduced and checked against Fig6_single_trial_robustness_reprocessed.csv (b and p, tolerance 1e-6) and the published
#         values (b, SE to 4 decimals, p to 3 decimals); N = 17 and trials entered / 1428 delivered (17 x 2 eye states x 42) = 93.0% / 92.0%. Any mismatch -> stop.
# Step 3: full model logMEP ~ P30_within + P30_between + resid_11_25_within + order_c + log_mu_within + (1 | subject).
# Step 4: VIF (performance::check_collinearity) and pooled Pearson correlations among P30_within, resid_11_25_within, order_c, log_mu_within.
# Output: TableS5_fullmodel_results.csv (this folder).
suppressMessages({library(lme4); library(lmerTest); library(dplyr); library(performance)})
MF <- paste0(data_dir, "/Revision_2026/Manuscript_Figures"); setwd(MF)
R8 <- file.path(MF, "MEP_reprocess_2026-10-22/main_qc2_008/analysis"); OUTF <- file.path(MF, "MEP_reprocess_2026-10-22/qc/v2/code008/TableS5_fullmodel_results.csv")
ctl <- lmerControl(optimizer = "bobyqa"); DELIVERED <- 17 * 2 * 42

X <- read.csv(file.path(R8, "Single_trial_robustness_input_reprocessed.csv"))
A <- X %>% filter(!is.na(mep_int_B), mep_int_B == intensity, mep_accept_B == 1) %>%
  mutate(present = as.integer(mep_present_B == 1), subject = factor(sub)) %>% group_by(id) %>% mutate(order_c = trial_order - mean(trial_order)) %>% ungroup()
wc <- function(d, v) { d %>% group_by(subject) %>% mutate("{v}_between" := mean(.data[[v]]), "{v}_within" := .data[[v]] - mean(.data[[v]])) %>% ungroup() }
prep <- function(d) { for (v in c("P30", "resid_11_25", "log_mu")) d <- wc(d, v); d }

MODELS <- list("(5) within/between, random intercept (main)" = "P30_within + P30_between", "(2) + residual 11-25 ms" = "P30_within + P30_between + resid_11_25_within",
               "(3) + trial order" = "P30_within + P30_between + order_c", "(6) + mu power" = "P30_within + P30_between + log_mu_within")
# published P30_within values (Table S5): b, SE, p (NA = reported only as < .001)
PUB <- data.frame(intensity = rep(c(100, 120), each = 4), model = rep(names(MODELS), 2),
                  b = c(0.0327, 0.0317, 0.0308, 0.0319, 0.0115, 0.0116, 0.0115, 0.0114), se = c(0.0069, 0.0070, 0.0069, 0.0070, 0.0045, 0.0045, 0.0045, 0.0045),
                  p = c(NA, NA, NA, NA, 0.011, 0.010, 0.011, 0.011), stringsAsFactors = FALSE)
PCT <- c("100" = 93.0, "120" = 92.0)
OLDD <- read.csv(file.path(R8, "Fig6_single_trial_robustness_reprocessed.csv"))

fit <- function(f, d) {   # fit and collect all warnings / convergence messages
  W <- character(0)
  m <- withCallingHandlers(lmer(f, d, REML = FALSE, control = ctl), warning = function(w) { W <<- c(W, conditionMessage(w)); invokeRestart("muffleWarning") })
  W <- c(W, m@optinfo$conv$lme4$messages, m@optinfo$warnings); if (isSingular(m)) W <- c(W, "singular fit")
  list(m = m, warn = if (length(W)) paste(unique(unlist(W)), collapse = " | ") else "none")
}
coefs <- function(m) { cf <- as.data.frame(summary(m)$coefficients); cf$term <- rownames(cf); cf }

MISM <- character(0); OUT <- data.frame()
for (lv in c(100, 120)) {
  u <- prep(A %>% filter(intensity == lv) %>% mutate(logMEP = log(mep_amp_B)))
  ns <- length(unique(u$subject)); pct <- 100 * nrow(u) / DELIVERED
  if (ns != 17) MISM <- c(MISM, sprintf("%d%%: %d participants (expected 17)", lv, ns))
  if (sprintf("%.1f", pct) != sprintf("%.1f", PCT[[as.character(lv)]])) MISM <- c(MISM, sprintf("%d%%: trials entered %.1f%% (expected %.1f%%)", lv, pct, PCT[[as.character(lv)]]))
  # ---- step 2: reproduction
  for (mn in names(MODELS)) {
    r <- fit(as.formula(paste("logMEP ~", MODELS[[mn]], "+ (1 | subject)")), u); cf <- coefs(r$m); cf <- cf[cf$term != "(Intercept)", ]
    o <- OLDD[OLDD$intensity == as.character(lv) & OLDD$model == mn, ]; oi <- match(cf$term, o$term)
    if (any(is.na(oi)) || any(abs(o$b[oi] - cf$Estimate) > 1e-6) || any(abs(o$p[oi] - cf$`Pr(>|t|)`) > 1e-6)) MISM <- c(MISM, sprintf("%d%% %s: differs from Fig6_single_trial_robustness_reprocessed.csv", lv, mn))
    pw <- cf[cf$term == "P30_within", ]; pb <- PUB[PUB$intensity == lv & PUB$model == mn, ]
    if (sprintf("%.4f", pw$Estimate) != sprintf("%.4f", pb$b)) MISM <- c(MISM, sprintf("%d%% %s: b %.6f vs published %.4f", lv, mn, pw$Estimate, pb$b))
    if (sprintf("%.4f", pw$`Std. Error`) != sprintf("%.4f", pb$se)) MISM <- c(MISM, sprintf("%d%% %s: SE %.6f vs published %.4f", lv, mn, pw$`Std. Error`, pb$se))
    if (is.na(pb$p) && pw$`Pr(>|t|)` >= .001) MISM <- c(MISM, sprintf("%d%% %s: p %.5f vs published < .001", lv, mn, pw$`Pr(>|t|)`))
    if (!is.na(pb$p) && sprintf("%.3f", pw$`Pr(>|t|)`) != sprintf("%.3f", pb$p)) MISM <- c(MISM, sprintf("%d%% %s: p %.5f vs published %.3f", lv, mn, pw$`Pr(>|t|)`, pb$p))
    OUT <- rbind(OUT, data.frame(section = "step2_reproduction", intensity = lv, model = mn, term = cf$term, var2 = NA, b = cf$Estimate, SE = cf$`Std. Error`, df = cf$df,
                                 t = cf$`t value`, p = cf$`Pr(>|t|)`, VIF = NA, r = NA, n_participants = ns, n_trials = nobs(r$m), pct_of_1428 = pct, warnings = r$warn))
  }
}
if (length(MISM)) { cat("STOP: reproduction failed\n", paste(MISM, collapse = "\n"), "\n"); stop("step 2 mismatch; steps 3-5 not run") }
cat("step 2: all reproduced (file values within 1e-6; published values to the reported decimals; N = 17; 93.0% / 92.0%)\n")

# ---- steps 3-4: full model, VIF, correlations
CV <- c("P30_within", "resid_11_25_within", "order_c", "log_mu_within")
for (lv in c(100, 120)) {
  u <- prep(A %>% filter(intensity == lv) %>% mutate(logMEP = log(mep_amp_B)))
  r <- fit(logMEP ~ P30_within + P30_between + resid_11_25_within + order_c + log_mu_within + (1 | subject), u); cf <- coefs(r$m)
  vf <- as.data.frame(check_collinearity(r$m)); vif <- setNames(vf$VIF, vf$Term)
  OUT <- rbind(OUT, data.frame(section = "step3_full_model", intensity = lv, model = "full", term = cf$term, var2 = NA, b = cf$Estimate, SE = cf$`Std. Error`, df = cf$df,
                               t = cf$`t value`, p = cf$`Pr(>|t|)`, VIF = unname(vif[cf$term]), r = NA, n_participants = length(unique(u$subject)), n_trials = nobs(r$m),
                               pct_of_1428 = 100 * nobs(r$m) / DELIVERED, warnings = r$warn))
  cm <- cor(u[, CV]); pr <- t(combn(CV, 2))
  OUT <- rbind(OUT, data.frame(section = "step4_correlation", intensity = lv, model = "pooled Pearson, within-participant variables", term = pr[, 1], var2 = pr[, 2], b = NA, SE = NA,
                               df = NA, t = NA, p = NA, VIF = NA, r = cm[pr], n_participants = length(unique(u$subject)), n_trials = nrow(u), pct_of_1428 = 100 * nrow(u) / DELIVERED, warnings = NA))
}
write.csv(OUT, OUTF, row.names = FALSE)

# ---- step 5: plain-text summary
fp <- function(p) if (p < .001) "< .001" else sub("^0", "", sprintf("%.3f", p))
for (lv in c(100, 120)) {
  FF <- OUT[OUT$section == "step3_full_model" & OUT$intensity == lv & OUT$term != "(Intercept)", ]
  cat(sprintf("\n%d%% RMT: N = %d participants, %d trials; warnings: %s\n", lv, FF$n_participants[1], FF$n_trials[1], FF$warnings[1]))
  for (i in seq_len(nrow(FF))) cat(sprintf("  %-20s b = %.4f, SE = %.4f, df = %.1f, t = %.3f, p %s, VIF = %.3f\n", FF$term[i], FF$b[i], FF$SE[i], FF$df[i], FF$t[i],
                                          ifelse(FF$p[i] < .001, "< .001", paste("=", fp(FF$p[i]))), FF$VIF[i]))
  Cr <- OUT[OUT$section == "step4_correlation" & OUT$intensity == lv, ]
  for (i in seq_len(nrow(Cr))) cat(sprintf("  r(%s, %s) = %.2f\n", Cr$term[i], Cr$var2[i], Cr$r[i]))
}
cat("\nwritten:", OUTF, "\n")
