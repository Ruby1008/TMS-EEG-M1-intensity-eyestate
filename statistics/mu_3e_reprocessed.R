data_dir <- "path/to/data"   # release: set this path
# qc/v2/make_code008.py copy (part 8): spTMS_008_Real_EC no longer excluded.
# copy made by MEP_reprocess_2026-10-22/code/make_reprocessed_copies.py (2026-10-22); MEPSET = main | excl_artifact; outputs in MEP_reprocess_2026-10-22/<MEPSET>/analysis/.
# section 3e of code/eyestate_null_round8.R (lines from '# ---- 3e mechanism' to the final bind_rows), unchanged except the MEP map
suppressMessages({library(lme4); library(lmerTest); library(dplyr)})
setwd(paste0(data_dir, "/Revision_2026/Manuscript_Figures"))
SET <- Sys.getenv("MEPSET", "main"); RD <- file.path(paste0(data_dir, "/Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22"), SET); OUTD <- file.path(RD, "analysis"); dir.create(OUTD, showWarnings = FALSE)

out <- data.frame(); ctl <- lmerControl(optimizer = "bobyqa")
add <- function(...) { r <- data.frame(..., stringsAsFactors = FALSE); out <<- bind_rows(out, r) }
# ---- 3e mechanism
pw <- read.csv("data/Prestim_power_single_trial.csv"); st <- read.csv("data/ERM_single_trial.csv"); mp <- read.csv(file.path(OUTD, "EEG_MEP_trial_map_reprocessed.csv"))
d <- merge(merge(pw, st[, c("id", "epoch", "P30")], by = c("id", "epoch")), mp[, c("id", "epoch", "sub", "eye", "mep_int_B", "mep_amp_B", "mep_accept_B")], by = c("id", "epoch"))
d$subject <- factor(d$sub)
for (lv in c(80, 100, 120)) for (xv in c("log_mu", "log_occ_alpha")) {
  xl <- ifelse(xv == "log_mu", "mu (FC3/FC1/C3/C1)", "occipital alpha (O1/O2/POz)")
  u <- d %>% filter(intensity == lv, !is.na(mep_int_B), mep_int_B == intensity, mep_accept_B == 1, mep_amp_B > 0) %>% mutate(logMEP = log(mep_amp_B)) %>%
    group_by(subject) %>% mutate(xw = .data[[xv]] - mean(.data[[xv]])) %>% ungroup()
  f <- lmer(logMEP ~ xw + (1 | subject), u, REML = FALSE, control = ctl); cf <- summary(f)$coefficients["xw", ]; sdx <- sd(u$xw)
  add(section = "3e", analysis = "single-trial log MEP ~ power_within + (1 | subject)", measure = xl, effect = paste0(lv, "% RMT"), n = length(unique(u$subject)), n_obs = nrow(u),
      b = cf["Estimate"], se = cf["Std. Error"], p = cf["Pr(>|t|)"], sd_x = sdx, change_per_SD = 100 * (exp(cf["Estimate"] * sdx) - 1), change_unit = "% MEP")
  v <- d %>% filter(intensity == lv) %>% group_by(subject) %>% mutate(xw = .data[[xv]] - mean(.data[[xv]])) %>% ungroup()
  g <- lmer(P30 ~ xw + (1 | subject), v, REML = FALSE, control = ctl); cg <- summary(g)$coefficients["xw", ]; sdv <- sd(v$xw)
  add(section = "3e", analysis = "single-trial P30 ~ power_within + (1 | subject)", measure = xl, effect = paste0(lv, "% RMT"), n = length(unique(v$subject)), n_obs = nrow(v),
      b = cg["Estimate"], se = cg["Std. Error"], p = cg["Pr(>|t|)"], sd_x = sdv, change_per_SD = cg["Estimate"] * sdv, change_unit = "uV P30")
}
write.csv(out, file.path(OUTD, "Mu_single_trial_3e_reprocessed.csv"), row.names = FALSE); print(out, digits = 3)
