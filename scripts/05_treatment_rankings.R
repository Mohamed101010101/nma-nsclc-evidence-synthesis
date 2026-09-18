# ==============================================================================
# Script: 05_treatment_rankings.R
# Purpose: Treatment Ranking Analysis via P-Scores & Hierarchy Visualization
# Outputs: outputs/tables/treatment_rankings.csv
#          outputs/figures/03_pscore_ranking.png (300 DPI Publication Figure)
# Package: netmeta & ggplot2
# ==============================================================================

suppressPackageStartupMessages({
  library(netmeta)
  library(ggplot2)
  library(scales)
})

cat("\n======================================================================\n")
cat(" [ANALYSIS 3/7] TREATMENT RANKING VIA P-SCORES & SUCRA ANALOGUE\n")
cat("======================================================================\n")

# 1. Load Clinical Trial Contrast Data
data_path <- "data/nsclc_trial_contrasts.csv"
if (!file.exists(data_path)) {
  stop(sprintf("Data file not found at: %s. Please run scripts/02_generate_data.R first.", data_path))
}
dat <- read.csv(data_path, stringsAsFactors = FALSE)
cat(sprintf(" - Loaded contrast dataset: %d comparisons across %d trials\n",
            nrow(dat), length(unique(dat$studlab))))

# 2. Fit Frequentist Graph-Theoretical Model
nma <- netmeta(
  TE = TE,
  seTE = seTE,
  treat1 = treat1,
  treat2 = treat2,
  studlab = studlab,
  data = dat,
  sm = "HR",
  reference.group = "Chemo",
  common = TRUE,
  random = TRUE,
  tol.multiarm = 0.005,
  details.chkmultiarm = FALSE
)

# 3. Compute P-Scores (Frequentist Analogue to Bayesian SUCRA)
# small.values = "good" because lower Hazard Ratio represents superior overall survival
rk <- netrank(nma, small.values = "good")
pscores_rand <- rk$ranking.random
trt_order <- names(sort(pscores_rand, decreasing = TRUE))

trt_labels_map <- c(
  "IO_Chemo"  = "IO + Chemo",
  "TKI_Chemo" = "TKI + Chemo",
  "Dual_IO"   = "Dual IO",
  "IO_Mono"   = "IO Monotherapy",
  "TKI"       = "TKI Monotherapy",
  "Chemo"     = "Chemotherapy"
)

# 4. Construct Comprehensive Treatment Ranking Table
df_rankings <- data.frame(
  Treatment = trt_labels_map[trt_order],
  Code = trt_order,
  Rank = 1:length(trt_order),
  Pscore_Random = round(pscores_rand[trt_order], 4),
  Pscore_Common = round(rk$ranking.common[trt_order], 4),
  HR_vs_Chemo_Random = ifelse(trt_order == "Chemo", "1.00 (Reference)",
                              sprintf("%.2f [%.2f; %.2f]", 
                                      exp(nma$TE.random[trt_order, "Chemo"]),
                                      exp(nma$lower.random[trt_order, "Chemo"]),
                                      exp(nma$upper.random[trt_order, "Chemo"]))),
  Pval_vs_Chemo = ifelse(trt_order == "Chemo", "Reference",
                         ifelse(nma$pval.random[trt_order, "Chemo"] < 0.0001, "< 0.0001",
                                sprintf("%.4f", nma$pval.random[trt_order, "Chemo"]))),
  stringsAsFactors = FALSE
)

# Export Ranking Table to CSV
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
output_tbl <- "outputs/tables/treatment_rankings.csv"
write.csv(df_rankings[, c("Treatment", "Rank", "Pscore_Random", "Pscore_Common", "HR_vs_Chemo_Random", "Pval_vs_Chemo")], 
          output_tbl, row.names = FALSE)
cat(sprintf(" - Exported ranking table to: %s\n", output_tbl))
print(df_rankings[, c("Treatment", "Rank", "Pscore_Random", "HR_vs_Chemo_Random", "Pval_vs_Chemo")])

# 5. Render Publication Ranking Hierarchy Bar Chart (300 DPI)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
output_fig <- "outputs/figures/03_pscore_ranking.png"
cat(sprintf("\n - Rendering Figure 3 to: %s ...\n", output_fig))

df_plot_rank <- df_rankings
df_plot_rank$Treatment <- factor(df_plot_rank$Treatment, levels = rev(df_rankings$Treatment))

colors_by_trt <- c(
  "Chemotherapy"    = "#718096",
  "IO Monotherapy"  = "#3182CE",
  "IO + Chemo"      = "#2B6CB0",
  "Dual IO"         = "#805AD5",
  "TKI Monotherapy" = "#D69E2E",
  "TKI + Chemo"     = "#DD6B20"
)

p_rank <- ggplot(df_plot_rank, aes(x = Pscore_Random, y = Treatment, fill = Treatment)) +
  geom_col(width = 0.65, alpha = 0.9, color = "#2D3748", linewidth = 0.4) +
  geom_text(aes(label = sprintf("Rank #%d | P-score: %.1f%%", Rank, Pscore_Random * 100)),
            hjust = -0.08, size = 4.2, fontface = "bold", color = "#1A202C") +
  scale_fill_manual(values = colors_by_trt) +
  scale_x_continuous(limits = c(0, 1.25), breaks = seq(0, 1, 0.2), 
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Treatment Ranking Hierarchy: Surface Under Cumulative Ranking (P-Scores)",
    subtitle = "Overall Survival in Advanced NSCLC (Frequentist random-effects model)",
    x = "P-Score (Certainty of Superiority over Competing Regimens)",
    y = NULL,
    caption = "P-score ranges from 0 (certain worst) to 1 (certain best).\nComputed using netrank(..., small.values = 'good')."
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "#E2E8F0", linetype = "dashed"),
    plot.title = element_text(face = "bold", size = 16, color = "#1A365D"),
    plot.subtitle = element_text(color = "#4A5568", size = 12, margin = margin(b = 15)),
    axis.text.y = element_text(face = "bold", size = 12, color = "#2D3748"),
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10))
  )

ggsave(output_fig, plot = p_rank, width = 10, height = 6.5, dpi = 300)

cat(sprintf(" [SUCCESS] Treatment rankings & Figure 3 saved cleanly: %s\n\n", output_fig))
