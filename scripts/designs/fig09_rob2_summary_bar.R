# ==============================================================================
# Design Script: scripts/designs/fig09_rob2_summary_bar.R
# Visual Target: Figure 9 - Cochrane RoB 2 Domain Summary Bar Chart
# Output File:   outputs/figures/09_rob2_summary_bar.png (300 DPI Publication Figure)
# Framework:     ggplot2 (Uses Cached RoB 2 Assessment Data)
# Standard:      Cochrane Guidelines & robvis Horizontal Stacked Percentage Plot
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

cat("\n======================================================================\n")
cat(" [DESIGN 9/10] FIGURE 9: COCHRANE RoB 2 DOMAIN SUMMARY BAR CHART\n")
cat("======================================================================\n")

# 1. Load RoB 2 Data
rob_rds  <- "outputs/models/rob2_data.rds"
data_csv <- "data/nsclc_rob2_assessments.csv"

if (!file.exists(rob_rds) || (file.exists(data_csv) && file.mtime(data_csv) > file.mtime(rob_rds))) {
  cat(" - RoB 2 cache missing or dataset modified. Running 06_rob2_analysis.R ...\n")
  source("scripts/analyses/06_rob2_analysis.R", local = new.env())
}

rob_bundle <- readRDS(rob_rds)
rob_df     <- rob_bundle$raw_data
n_studies  <- nrow(rob_df)

domains <- c("D1", "D2", "D3", "D4", "D5", "Overall")
domain_labels <- c(
  "D1"      = "D1: Randomisation process",
  "D2"      = "D2: Deviations from intended interventions",
  "D3"      = "D3: Missing outcome data",
  "D4"      = "D4: Measurement of the outcome",
  "D5"      = "D5: Selection of the reported result",
  "Overall" = "Overall Risk of Bias"
)

# 2. Reshape Frequencies for Stacked Horizontal Percentage Bar
levels_rob <- c("Low risk", "Some concerns", "High risk")
colors_rob <- c(
  "Low risk"      = "#38A169", # Emerald Green
  "Some concerns" = "#ECC94B", # Yellow / Amber
  "High risk"     = "#E53E3E"  # Crimson Red
)

bar_rows <- list()
for (d in domains) {
  counts <- table(factor(rob_df[[d]], levels = levels_rob))
  for (lev in levels_rob) {
    cnt <- as.numeric(counts[lev])
    pct <- cnt / n_studies * 100
    bar_rows[[length(bar_rows) + 1]] <- data.frame(
      Domain = d,
      Domain_Label = domain_labels[d],
      Judgment = lev,
      Count = cnt,
      Percentage = pct,
      stringsAsFactors = FALSE
    )
  }
}

df_bar <- do.call(rbind, bar_rows)
# Factor ordering: Overall at the bottom, D1 at the top (reversed for ggplot coord_flip / horizontal)
df_bar$Domain_Label <- factor(df_bar$Domain_Label, levels = rev(domain_labels))
df_bar$Judgment <- factor(df_bar$Judgment, levels = rev(levels_rob)) # Reversed so Low risk stacks first from left

# 3. Render Publication RoB 2 Summary Bar Chart (300 DPI)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
output_fig <- "outputs/figures/09_rob2_summary_bar.png"
cat(sprintf(" - Rendering Figure 9 to: %s ...\n", output_fig))

p_bar <- ggplot(df_bar, aes(x = Domain_Label, y = Percentage, fill = Judgment)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.65, 
           color = "#2D3748", linewidth = 0.35, alpha = 0.95) +
  geom_text(
    data = subset(df_bar, Percentage >= 6.0),
    aes(label = sprintf("%.1f%%", Percentage)),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    color = ifelse(subset(df_bar, Percentage >= 6.0)$Judgment == "Some concerns", "#1A202C", "#FFFFFF"),
    fontface = "bold",
    size = 3.8
  ) +
  coord_flip() +
  scale_fill_manual(
    name = "Risk of Bias Judgment:",
    values = colors_rob,
    breaks = levels_rob,
    labels = c("Low risk", "Some concerns", "High risk")
  ) +
  scale_y_continuous(
    limits = c(0, 100.1),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = c(0.01, 0.01)
  ) +
  labs(
    title = "Risk of Bias 2 (RoB 2) Domain-Level Summary Distribution",
    subtitle = sprintf("Proportional evaluation across %d Randomized Controlled Trials in First-Line NSCLC", n_studies),
    x = NULL,
    y = "Proportion of Clinical Trials (%)",
    caption = "Evaluated according to the revised Cochrane Risk of Bias tool for randomized trials (RoB 2)."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, color = "#1A365D", margin = margin(b = 6)),
    plot.subtitle = element_text(color = "#4A5568", size = 11, margin = margin(b = 14)),
    plot.caption = element_text(color = "#718096", size = 9, hjust = 0, margin = margin(t = 12)),
    axis.text.y = element_text(face = "bold", size = 11, color = "#2D3748"),
    axis.text.x = element_text(size = 10.5, color = "#4A5568"),
    axis.title.x = element_text(face = "bold", size = 11, color = "#1A365D", margin = margin(t = 8)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#E2E8F0", linetype = "dashed"),
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_text(face = "bold", size = 10.5, color = "#1A365D"),
    legend.text = element_text(size = 10, color = "#2D3748"),
    legend.background = element_rect(fill = "#F7FAFC", color = "#E2E8F0", linewidth = 0.5),
    legend.margin = margin(t = 4, b = 4, l = 10, r = 10),
    plot.margin = margin(16, 22, 14, 16)
  )

ggsave(output_fig, plot = p_bar, width = 10.5, height = 6.2, dpi = 300)

cat(sprintf(" [SUCCESS] Figure 9 (RoB 2 Summary Bar Chart) rendered cleanly: %s\n\n", output_fig))
