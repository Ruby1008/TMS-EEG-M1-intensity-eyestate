data_dir <- "path/to/data"   # release: set this path
# Round 15 merge (2026-10-22): ITPC rows of eyestate_null_round8.R 3c (ii) and (v) with ITPC rule v2 (input data/EyeState_allavailable_long_ITPC_v2.csv);
# models copied unchanged; appends to data/EyeState_null_robustness_ITPC_v2.csv (written by eyestate_itpc_v2.py).
suppressMessages({library(lme4); library(lmerTest); library(dplyr)})
setwd(paste0(data_dir, "/Revision_2026/Manuscript_Figures"))
R0 <- read.csv("data/EyeState_null_robustness_ITPC_v2.csv", check.names = FALSE)
L <- read.csv("data/EyeState_allavailable_long_ITPC_v2.csv") %>% mutate(subject = factor(sub), state = factor(eye, levels = c("EO", "EC")), intens = factor(intensity))
ORD <- read.csv("data/EyeState_block_order.csv"); L <- L %>% left_join(ORD[, c("sub", "order")], by = "sub")
contrasts(L$state) <- contr.sum(2); contrasts(L$intens) <- contr.sum(3)
out <- data.frame(); ctl <- lmerControl(optimizer = "bobyqa")
add <- function(...) { r <- data.frame(..., stringsAsFactors = FALSE); out <<- bind_rows(out, r) }
for (m in unique(L$measure)) {
  d <- L %>% filter(measure == m)
  f <- lmer(value ~ state * intens + (1 | subject), d, REML = FALSE, control = ctl); a <- as.data.frame(anova(f, type = 3))
  for (tm in c("state", "state:intens")) add(section = "3c", analysis = "(ii) LMM, all available participants", measure = m, effect = tm, n = length(unique(d$subject)), n_obs = nrow(d),
                                           F = a[tm, "F value"], df1 = a[tm, "NumDF"], df2 = a[tm, "DenDF"], p = a[tm, "Pr(>F)"])
  dm <- d %>% filter(main17 == "True" | main17 == TRUE)
  dm$order <- factor(dm$order); contrasts(dm$order) <- contr.sum(2)
  f2 <- lmer(value ~ state * intens * order + (1 | subject), dm, REML = FALSE, control = ctl); a2 <- as.data.frame(anova(f2, type = 3))
  for (tm in c("state", "state:order", "order")) add(section = "3c", analysis = "(v) block order as between-participant factor (LMM, main set)", measure = m, effect = tm, n = length(unique(dm$subject)),
                                                    n_EO_first = length(unique(dm$sub[dm$order == "EO first"])), F = a2[tm, "F value"], df1 = a2[tm, "NumDF"], df2 = a2[tm, "DenDF"], p = a2[tm, "Pr(>F)"])
}
R1 <- bind_rows(R0, out); write.csv(R1, "data/EyeState_null_robustness_ITPC_v2.csv", row.names = FALSE)
