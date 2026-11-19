data_dir <- "path/to/data"   # release: set this path
# _1114 copy (user 2026-11-14, Figures_2026-11-14_minor) of Figures_2026-11-08t/code/Fig2_1108t.R (run with MEPSET=main_qc2, T_MODE=rev).
# M1114=repro: unchanged -> _check/repro; M1114=rev: y axis 'MEP amplitude (mV)'. Lines marked # 1114.
# _1108t copy (user 2026-10-28, Figures_2026-11-08t) of MEP_reprocess_2026-10-22/qc/v2/Fig2_p6_qc2.R. Changes (T_MODE=rev): title two lines without semicolon, x-axis (and unused legend) labels '80% RMT (EC)' etc.; output Figure2_1108t.png/.pdf
# Part 5/6 copy of qc/figures_p6/Fig2_p6.R with output in qc/v2/figures (MEPSET=main_qc2); nothing else changed.
# Part 6 copy made by qc/figures_p6/make_p6_copies.py (2026-10-23): brackets only (see that file); MEPSET = main_qc for the 10-30/10-31 manuscript.
# A3 copy made by qc/figures_a3/make_a3_figure_copies.py (part 4, 2026-10-22).
# copy made by MEP_reprocess_2026-10-22/code/make_reprocessed_copies.py (2026-10-22); MEPSET = main | excl_artifact; outputs in MEP_reprocess_2026-10-22/<MEPSET>/analysis/.
# Fig2_r13 (2026-10-25, round 13 A-1): copy of Fig2_r9_MEP_boxplot.R; only change: title on two lines so that it is not cut off; output Fig2_r13.png/.tiff.
# Fig2_r9 (2026-10-21, round 9: title per user; otherwise identical to Fig3_new_MEP_boxplot.R). Fig3_new (2026-10-14): copy of MEP_Amplitude_EOEC_All23_Boxplot_20250605.R (from <data_dir>/spTMS_Formal_EMG, MEP.rar).
# Changes ONLY: (1) input = data/MEP_matrix_new.csv (trial_accept + non-delivered trials excluded; sub 8 EC_120 missing);
# (2) working dir / output file names; (3) font_import() and windows() removed (non-interactive run; fonts already imported);
# (4) 2026-10-19: complete cases only (n = 22, same participants in EO and EC; user decision 1-c);
# (7) 2026-10-20: y axis log10 (coord_trans; box statistics computed on mV), ticks 0.01/0.1/1 mV; stars = Holm paired t on log MEP;
# (5) red post hoc brackets added as in the published Fig 3 (drawn outside R in the original); stars from data/MEP_posthoc.csv
#     (Bonferroni Wilcoxon, eye state collapsed): *** p < .001, ** p < .01, * p < .05; (6) saved at 600 dpi (high resolution) only.
# Clear all
rm(list=ls())
SET <- Sys.getenv("MEPSET", "main"); RD <- file.path(paste0(data_dir, "/Revision_2026/Manuscript_Figures/MEP_reprocess_2026-10-22"), SET); OUTD <- file.path(RD, "analysis"); dir.create(OUTD, showWarnings = FALSE)


# 載入 extrafont 套件
library(extrafont)
TMODE <- Sys.getenv("T_MODE", "rev")   # 1108t: repro | rev
M1114 <- Sys.getenv("M1114", "rev"); stopifnot(M1114 %in% c("repro", "rev"), TMODE == "rev")   # 1114
# 匯入系統字型
# font_import(prompt = FALSE)  # removed: fonts already imported on this machine
loadfonts(device = "win")    # 適用 Windows


library(reshape2)
library(ggplot2)
library(dplyr)

setwd(paste0(data_dir, "/Revision_2026/Manuscript_Figures"))

# 讀取 CSV 檔案
data <- read.csv(file.path(RD, "MEP_matrix_reprocessed.csv"), header = TRUE); names(data)[names(data) == "subject"] <- "sub"
# 2026-10-19 (user decision 1-c): complete cases only (22 participants with all six cells), so that EO and EC have the same participants as the ANOVA
data <- data[complete.cases(data), ]; cat("n =", nrow(data), "\n")

# 資料轉 long format
data_long <- melt(
  data,
  id.vars = "sub",
  variable.name = "Condition",
  value.name = "Amplitude"
)

# 設定條件順序
data_long$Condition <- factor(
  data_long$Condition,
  levels = c("EC_80", "EO_80", "EC_100", "EO_100", "EC_120", "EO_120")
)

# 指定每個條件對應的 x 軸位置
# 調整為 1,2 | 3,4 | 5,6
# 維持每組 EC–EO 的間距 = 1
# 且縮短 80%(EO)→100%(EC)、100%(EO)→120%(EC) 的跨強度間距（由 2 縮為 1）
lookup <- c(
  "EC_80"  = 1, "EO_80"  = 2,
  "EC_100" = 3, "EO_100" = 4,
  "EC_120" = 5, "EO_120" = 6
)
data_long$xpos <- lookup[as.character(data_long$Condition)]
data_long <- data_long[!is.na(data_long$Amplitude), ]
print(table(data_long$Condition))

# 外部視窗（Windows）
# windows(width = 13, height = 10, pointsize = 22)  # removed: non-interactive
library(grid)  # for unit()

p <- ggplot(data_long, aes(x = xpos, y = Amplitude, fill = Condition)) +
  geom_boxplot(
    aes(group = Condition),
    color = "black",
    size = 1.1,         # box 線條加粗
    outlier.shape = 16,
    outlier.colour = "black",
    outlier.size = 2.8, # 離群點放大
    width = 0.7
  ) +
  scale_x_continuous(
    breaks = 1:6,
    labels = if (TMODE == "repro") c("80%RMT(EC)", "80%RMT(EO)", "100%RMT(EC)", "100%RMT(EO)", "120%RMT(EC)", "120%RMT(EO)") else
             c("80% RMT (EC)", "80% RMT (EO)", "100% RMT (EC)", "100% RMT (EO)", "120% RMT (EC)", "120% RMT (EO)")   # 1108t
  ) +
  scale_y_continuous(
    breaks = c(0.01, 0.1, 1), labels = c("0.01", "0.1", "1")   # 2026-10-20 (user decision 1-3): log axis via coord_trans below; box statistics on mV values
  ) +
  scale_fill_manual(
    values = c(
      "EO_80"  = rgb(249/255, 149/255,  78/255),
      "EO_100" = rgb(249/255, 114/255,  19/255),
      "EO_120" = rgb(197/255,  85/255,   5/255),
      "EC_80"  = rgb(122/255, 184/255, 184/255),
      "EC_100" = rgb(109/255, 173/255, 157/255),
      "EC_120" = rgb( 57/255, 106/255, 106/255)
    ),
    labels = if (TMODE == "repro") c("EC_80" = "80%RMT(EC)", "EO_80" = "80%RMT(EO)", "EC_100" = "100%RMT(EC)", "EO_100" = "100%RMT(EO)", "EC_120" = "120%RMT(EC)", "EO_120" = "120%RMT(EO)") else
             c("EC_80" = "80% RMT (EC)", "EO_80" = "80% RMT (EO)", "EC_100" = "100% RMT (EC)", "EO_100" = "100% RMT (EO)", "EC_120" = "120% RMT (EC)", "EO_120" = "120% RMT (EO)")   # 1108t (legend not shown)
  ) +
  labs(
    x = "Condition",
    y = if (M1114 == "repro") "MEP Amplitude (mV)" else "MEP amplitude (mV)",   # 1114
    fill = "Condition"
  ) +
  theme_minimal(base_family = "Times New Roman", base_size = 26) +  # 全域字級偏大
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "none",  # 關閉 legend
    panel.grid = element_blank(),
    axis.line = element_line(size = 1.2, colour = "black"),
    axis.ticks = element_line(size = 1.2, colour = "black"),
    axis.ticks.length = unit(7, "pt"),
    axis.text.x = element_text(
      size = 18, face = "bold", colour = "black",
      angle = 35, vjust = 1, hjust = 1           # x 軸刻度較小且傾斜
    ),
    axis.text.y = element_text(size = 22, face = "bold", colour = "black"),
    axis.title.x = element_text(size = 20, face = "bold", colour = "black", margin = margin(t = 16)),
    axis.title.y = element_text(size = 22, face = "bold", colour = "black", margin = margin(r = 16)),
    plot.margin = margin(20, 20, 30, 20)  # 多留底部空間避免刻度被切到
  )


# ⬇️ 高畫質輸出（600 dpi PNG；以目前工作目錄儲存）

# ---- red post hoc brackets (layout of the published Fig 3) --------------------------------------
# 2026-10-20 (user decision 1): stars from the main post hoc tests = paired t-tests on log MEP (eye states averaged), Holm (data/MEP_stats_A_vs_B.csv, version B)
ph <- read.csv(file.path(OUTD, "MEP_stats_A_vs_B_reprocessed.csv"), fileEncoding = "UTF-8-BOM"); ph <- ph[ph$version == "B" & grepl("post hoc", ph$block), ]
stars <- function(t) { t <- sub("Wilcoxon ", "", t); q <- ph$p_corr[ph$item == t]; if (q < .001) "***" else if (q < .01) "**" else if (q < .05) "*" else "n.s." }
red <- "#C00000"; ymax <- max(data_long$Amplitude, na.rm = TRUE); ymin <- min(data_long$Amplitude, na.rm = TRUE)
# log axis: bracket heights as multiplicative steps
y1 <- ymax * 10^0.15; y2 <- ymax * 10^0.42; y3 <- ymax * 10^0.69; y4 <- ymax * 10^0.96; legf <- 10^0.08   # part 6: three heights (y2, y3, y4)
leg <- 0   # unused on the log axis (bracket legs drawn multiplicatively)
br <- function(x0, x1, y) list(annotate("segment", x = x0, xend = x1, y = y, yend = y, colour = red, linewidth = 1.3),
                               annotate("segment", x = x0, xend = x0, y = y, yend = y / legf, colour = red, linewidth = 1.3),
                               annotate("segment", x = x1, xend = x1, y = y, yend = y / legf, colour = red, linewidth = 1.3))
# part 6: group centres (eye states averaged) 80% = 1.5, 100% = 3.5, 120% = 5.5; each bracket spans only its two groups, stars at the centre
p <- p + br(1.5, 3.5, y2) + br(3.5, 5.5, y3) + br(1.5, 5.5, y4) +
  annotate("text", x = 2.5, y = y2 * 10^0.07, label = stars("Wilcoxon 100 vs 80"), colour = red, size = 9, family = "Times New Roman") +
  annotate("text", x = 4.5, y = y3 * 10^0.07, label = stars("Wilcoxon 120 vs 100"), colour = red, size = 9, family = "Times New Roman") +
  annotate("text", x = 3.5, y = y4 * 10^0.07, label = stars("Wilcoxon 120 vs 80"), colour = red, size = 9, family = "Times New Roman") +
  coord_trans(y = "log10", ylim = c(ymin, ymax * 1.05), clip = "off") +
  theme(plot.margin = margin(40, 20, 30, 20)) +
  # 2026-10-21 (round 8): one-line conclusion title above the post hoc brackets
  labs(title = if (TMODE == "repro") "MEP amplitude increases with stimulation intensity;\nno significant effect of eye state" else
                "MEP amplitude increases with stimulation intensity\nwithout a significant effect of eye state") +   # 1108t
  theme(plot.title = element_text(family = "Times New Roman", face = "bold", size = 24, hjust = 0.5, margin = margin(b = 175)))   # part 6: room for the third bracket
A3D <- paste0(data_dir, "/Revision_2026/Manuscript_Figures/Figures_2026-11-14_minor")   # 1114: output folder
STEM <- if (M1114 == "repro") file.path(A3D, "_check", "repro", "repro_Figure2_1108t") else file.path(A3D, "Figure2_1114")   # 1114
if (TMODE == "repro") {
  ggsave(file.path(A3D, "_check", "_scratch", "repro_Fig2_p6.png"), plot = p, width = 13, height = 11.5, units = "in", dpi = 600)
} else {
  ggsave(paste0(STEM, ".png"), plot = p, width = 13, height = 11.5, units = "in", dpi = 600)   # 1114
  ggsave(paste0(STEM, ".pdf"), plot = p, width = 13, height = 11.5, units = "in", device = cairo_pdf)
}
