# ==============================================================================
# Design Script: scripts/designs/fig08_rob2_traffic_light.R
# Visual Target: Figure 8 - Cochrane RoB 2 Study-Level Traffic Light Matrix
# Output File:   outputs/figures/08_rob2_traffic_light.png (300 DPI Publication Figure)
# Framework:     ggplot2 (Uses Cached RoB 2 Assessment Data)
# Standard:      Cochrane RoB 2 Guidelines & robvis Publication Aesthetics
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
})

cat("\n======================================================================\n")
cat(" [DESIGN 8/10] FIGURE 8: COCHRANE RoB 2 TRAFFIC LIGHT MATRIX PLOT\n")
cat("======================================================================\n")

# 1. Load RoB 2 Data
rob_rds  <- "outputs/models/rob2_data.rds"
data_csv <- "data/nsclc_rob2_assessments.csv"

if (!file.exists(rob_rds) || (file.exists(data_csv) && file.mtime(data_csv) > file.mtime(rob_rds))) {
  cat(" - RoB 2 cache missing or dataset modified. Running 06_rob2_analysis.R ...\n")
  source("scripts/analyses/06_rob2_analysis.R", local = new.env())
}

rob_bundle <- readRDS(rob_rds)
rob_df <- rob_bundle$raw_data
cat(sprintf(" - Loaded RoB 2 assessments for %d trials.\n", nrow(rob_df)))

# 2. Reshape into Long Format for High-Precision ggplot2 Matrix
domains <- c("D1", "D2", "D3", "D4", "D5", "Overall")
domain_labels <- c(
  "D1"      = "D1: Randomisation\nprocess",
  "D2"      = "D2: Deviations from\nintended intervention",
  "D3"      = "D3: Missing\noutcome data",
  "D4"      = "D4: Measurement of\nthe outcome",
  "D5"      = "D5: Selection of\nreported result",
  "Overall" = "Overall\nRisk of Bias"
)

# Sort studies by year then name (reverse order for top-to-bottom reading in ggplot2)
rob_df_sorted <- rob_df[order(rob_df$year, rob_df$studlab, decreasing = TRUE), ]
study_levels <- rob_df_sorted$studlab

plot_rows <- list()
for (i in 1:nrow(rob_df_sorted)) {
  st <- rob_df_sorted$studlab[i]
  for (d in domains) {
    val <- rob_df_sorted[[d]][i]
    symbol <- switch(val,
      "Low risk"      = "+",
      "Some concerns" = "?",
      "High risk"     = "-",
      "?"
    )
    plot_rows[[length(plot_rows) + 1]] <- data.frame(
      Study = st,
      Domain = d,
      Domain_Label = domain_labels[d],
      Judgment = val,
      Symbol = symbol,
      Is_Overall = (d == "Overall"),
      Row_Index = i,
      stringsAsFactors = FALSE
    )
  }
}

df_plot <- do.call(rbind, plot_rows)
df_plot$Study <- factor(df_plot$Study, levels = rev(study_levels))
df_plot$Domain_Label <- factor(df_plot$Domain_Label, levels = domain_labels)
df_plot$Judgment <- factor(df_plot$Judgment, levels = c("Low risk", "Some concerns", "High risk"))

# 3. Create Alternating Zebra Row Shading
zebra_df <- data.frame(
  Study = rev(study_levels),
  Y_Index = seq_along(study_levels),
  Fill = ifelse(seq_along(study_levels) %% 2 == 0, "#F7FAFC", "#FFFFFF"),
  stringsAsFactors = FALSE
)

# 4. Color Palette & Aesthetic Configuration
colors_rob <- c(
  "Low risk"      = "#38A169", # Classic Cochrane Emerald Green
  "Some concerns" = "#ECC94B", # Cochrane Yellow / Amber
  "High risk"     = "#E53E3E"  # Cochrane Crimson Red
)

symbol_colors <- c(
  "Low risk"      = "#FFFFFF",
  "Some concerns" = "#2D3748",
  "High risk"     = "#FFFFFF"
)

# 5. Render Publication Traffic Light Matrix (300 DPI)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
output_fig <- "outputs/figures/08_rob2_traffic_light.png"
cat(sprintf(" - Rendering Figure 8 to: %s ...\n", output_fig))

p_tl <- ggplot(df_plot, aes(x = Domain_Label, y = Study)) +
  # Zebra row striping
  geom_tile(aes(fill = NULL), width = 1, height = 1, fill = ifelse(as.numeric(df_plot$Study) %% 2 == 0, "#F8FAFC", "#FFFFFF")) +
  # Vertical divider separating D5 and Overall
  geom_vline(xintercept = 5.5, color = "#CBD5E0", linewidth = 1.1, linetype = "solid") +
  # Grid lines
  geom_hline(yintercept = seq(0.5, length(study_levels) + 0.5, 1), color = "#EDF2F7", linewidth = 0.5) +
  geom_vline(xintercept = seq(0.5, 6.5, 1), color = "#EDF2F7", linewidth = 0.5) +
  # Circular Traffic Light Glyphs
  geom_point(aes(color = Judgment), size = 9.5) +
  scale_color_manual(
    name = "Risk of Bias Judgment:",
    values = colors_rob,
    labels = c("Low risk (+)", "Some concerns (?)", "High risk (-)")
  ) +
  # Overlay contrasting symbol text
  geom_text(data = subset(df_plot, Judgment == "Low risk"),
            aes(label = "+"), color = "#FFFFFF", fontface = "bold", size = 4.8) +
  geom_text(data = subset(df_plot, Judgment == "Some concerns"),
            aes(label = "?"), color = "#1A202C", fontface = "bold", size = 4.2) +
  geom_text(data = subset(df_plot, Judgment == "High risk"),
            aes(label = "-"), color = "#FFFFFF", fontface = "bold", size = 5.0) +
  labs(
    title = "Cochrane Risk of Bias 2 (RoB 2) Traffic Light Evaluation",
    subtitle = "Trial-Level Methodological Quality Across 24 Randomized Controlled Trials in First-Line NSCLC",
    x = NULL,
    y = NULL,
    caption = paste0(
      "Domains: D1 = Bias arising from the randomization process; D2 = Bias due to deviations from intended interventions;\n",
      "D3 = Bias due to missing outcome data; D4 = Bias in measurement of the outcome; D5 = Bias in selection of the reported result.\n",
      "Visualized in accordance with Cochrane RoB 2 criteria (Sterne et al. BMJ 2019)."
    )
  ) +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16, color = "#1A365D", margin = margin(b = 6)),
    plot.subtitle = element_text(color = "#4A5568", size = 11.5, margin = margin(b = 16)),
    plot.caption = element_text(color = "#718096", size = 9.5, hjust = 0, lineheight = 1.3, margin = margin(t = 14)),
    axis.text.x.top = element_text(face = "bold", size = 11, color = "#2D3748", lineheight = 1.1),
    axis.text.y = element_text(face = "bold", size = 10.5, color = "#1A202C"),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11, color = "#1A365D"),
    legend.text = element_text(size = 11, color = "#2D3748"),
    legend.background = element_rect(fill = "#F7FAFC", color = "#E2E8F0", linewidth = 0.5),
    legend.margin = margin(t = 6, b = 6, l = 12, r = 12),
    plot.margin = margin(20, 24, 16, 24)
  )

ggsave(output_fig, plot = p_tl, width = 11.5, height = 12.5, dpi = 300)

cat(sprintf(" [SUCCESS] Figure 8 (RoB 2 Traffic Light Matrix) rendered cleanly: %s\n\n", output_fig))
