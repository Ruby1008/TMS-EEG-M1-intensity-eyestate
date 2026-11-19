data_dir <- "path/to/data"   # release: set this path
# Round 8 analyses (user 2026-10-20 evening, section 1), part B (R); appends to data/EyeState_null_robustness.csv written by eyestate_null_round8.py.
# 3c (ii) LMM with all available participants: value ~ state * intensity + (1 | subject) for log MEP (all MEP participants), ERM P30/N45/P60/N100 and IMF6/IMF7 ITPC
#      (longest window) (EO 20, EC 18 recordings); ML, sum contrasts, Type III F (Satterthwaite).
# 3c (v) block order (EO first / EC first, from the recording start times) as a between-participant factor: value ~ state * intensity * order + (1 | subject) on the main
#      sets (17 EEG participants; 22 MEP participants); state x order reported.
# 3e mechanism (single trial, per intensity 80/100/120% RMT; 17 participants; trials as Figure 6 for MEP, all epochs for P30):
#      log MEP ~ X_within + (1 | subject) and P30 ~ X_within + (1 | subject), X = pre-stimulus mu (FC3/FC1/C3/C1) or occipital alpha (O1/O2/POz) log10 8-12 Hz power,
#      centred within participant (per intensity); b, p and the change per 1 SD of X_within (MEP: %, P30: uV).
suppressMessages({library(lme4); library(lmerTest); library(dplyr)})
setwd(paste0(data_dir, "/Revision_2026/Manuscript_Figures"))
R0 <- read.csv("data/EyeState_null_robustness.csv", check.names = FALSE)
L <- read.csv("data/EyeState_allavailable_long.csv") %>% mutate(subject = factor(sub), state = factor(eye, levels = c("EO", "EC")), intens = factor(intensity))
ORD <- read.csv("data/EyeState_block_order.csv"); L <- L %>% left_join(ORD[, c("sub", "order")], by = "sub")
contrasts(L$state) <- contr.sum(2); contrasts(L$intens) <- contr.sum(3)
out <- data.frame(); ctl <- lmerControl(optimizer = "bobyqa")
add <- function(...) { r <- data.frame(..., stringsAsFactors = FALSE); out <<- bind_rows(out, r) }
for (m in unique(L$measure)) {
  d <- L %>% filter(measure == m)
  f <- lmer(value ~ state * intens + (1 | subject), d, REML = FALSE, control = ctl); a <- as.data.frame(anova(f, type = 3))
  for (tm in c("state", "state:intens")) add(section = "3c", analysis = "(ii) LMM, all available participants", measure = m, effect = tm, n = length(unique(d$subject)), n_obs = nrow(d),
                                           F = a[tm, "F value"], df1 = a[tm, "NumDF"], df2 = a[tm, "DenDF"], p = a[tm, "Pr(>F)"])
  dm <- if (m == "log MEP") d %>% group_by(sub) %>% filter(n() == 6) %>% ungroup() else d %>% filter(main17 == "True" | main17 == TRUE)
  dm$order <- factor(dm$order); contrasts(dm$order) <- contr.sum(2)
  f2 <- lmer(value ~ state * intens * order + (1 | subject), dm, REML = FALSE, control = ctl); a2 <- as.data.frame(anova(f2, type = 3))
  for (tm in c("state", "state:order", "order")) add(section = "3c", analysis = "(v) block order as between-participant factor (LMM, main set)", measure = m, effect = tm, n = length(unique(dm$subject)),
                                                    n_EO_first = length(unique(dm$sub[dm$order == "EO first"])), F = a2[tm, "F value"], df1 = a2[tm, "NumDF"], df2 = a2[tm, "DenDF"], p = a2[tm, "Pr(>F)"])
}
# ---- 3e mechanism
pw <- read.csv("data/Prestim_power_single_trial.csv"); st <- read.csv("data/ERM_single_trial.csv"); mp <- read.csv("data/EEG_MEP_trial_map.csv")
d <- merge(merge(pw, st[, c("id", "epoch", "P30")], by = c("id", "epoch")), mp[, c("id", "epoch", "sub", "eye", "mep_int_B", "mep_amp_B", "mep_accept_B")], by = c("id", "epoch"))
d$subject <- factor(d$sub)
for (lv in c(80, 100, 120)) for (xv in c("log_mu", "log_occ_alpha")) {
  xl <- ifelse(xv == "log_mu", "mu (FC3/FC1/C3/C1)", "occipital alpha (O1/O2/POz)")
  u <- d %>% filter(intensity == lv, id != "spTMS_008_Real_EC", !is.na(mep_int_B), mep_int_B == intensity, mep_accept_B == 1, mep_amp_B > 0) %>% mutate(logMEP = log(mep_amp_B)) %>%
    group_by(subject) %>% mutate(xw = .data[[xv]] - mean(.data[[xv]])) %>% ungroup()
  f <- lmer(logMEP ~ xw + (1 | subject), u, REML = FALSE, control = ctl); cf <- summary(f)$coefficients["xw", ]; sdx <- sd(u$xw)
  add(section = "3e", analysis = "single-trial log MEP ~ power_within + (1 | subject)", measure = xl, effect = paste0(lv, "% RMT"), n = length(unique(u$subject)), n_obs = nrow(u),
      b = cf["Estimate"], se = cf["Std. Error"], p = cf["Pr(>|t|)"], sd_x = sdx, change_per_SD = 100 * (exp(cf["Estimate"] * sdx) - 1), change_unit = "% MEP")
  v <- d %>% filter(intensity == lv) %>% group_by(subject) %>% mutate(xw = .data[[xv]] - mean(.data[[xv]])) %>% ungroup()
  g <- lmer(P30 ~ xw + (1 | subject), v, REML = FALSE, control = ctl); cg <- summary(g)$coefficients["xw", ]; sdv <- sd(v$xw)
  add(section = "3e", analysis = "single-trial P30 ~ power_within + (1 | subject)", measure = xl, effect = paste0(lv, "% RMT"), n = length(unique(v$subject)), n_obs = nrow(v),
      b = cg["Estimate"], se = cg["Std. Error"], p = cg["Pr(>|t|)"], sd_x = sdv, change_per_SD = cg["Estimate"] * sdv, change_unit = "uV P30")
}
R1 <- bind_rows(R0 %>% mutate(across(everything(), as.character)), out %>% mutate(across(everything(), as.character)))
write.csv(R1, "data/EyeState_null_robustness.csv", row.names = FALSE); print(out, digits = 3)
