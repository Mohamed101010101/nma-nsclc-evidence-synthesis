# ==============================================================================
# Design Script: scripts/designs/fig10_cinema_traffic_light.R
# Visual Target: Figure 10 - CINeMA Confidence Traffic Light & Evidence Matrix
# Output File:   outputs/figures/10_cinema_traffic_light.png (300 DPI Publication Figure)
# Framework:     ggplot2 & CINeMA Framework (Salanti et al., Nikolakopoulou et al. 2020)
# Evaluation:    6 Methodological Domains + Overall Confidence across 15 Comparisons
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
})

cat("\n======================================================================\n")
cat(" [DESIGN 10/10] FIGURE 10: CINeMA CONFIDENCE TRAFFIC LIGHT MATRIX\n")
cat("======================================================================\n")

# 1. Load CINeMA Evaluation Data
cinema_rds <- "outputs/models/cinema_data.rds"
cinema_csv <- "outputs/tables/cinema_summary_table.csv"

if (!file.exists(cinema_rds) || (file.exists(cinema_csv) && file.mtime(cinema_csv) > file.mtime(cinema_rds))) {
  cat(" - CINeMA cache missing or table modified. Running 07_cinema_evaluation.R ...\n")
  source("scripts/analyses/07_cinema_evaluation.R", local = new.env())
}

cinema_obj <- readRDS(cinema_rds)
df_cin <- cinema_obj$data
cat(sprintf(" - Loaded CINeMA data: %d pairwise comparisons.\n", nrow(df_cin)))

# 2. Reshape into Grid for Precise Visualization
# Columns in matrix:
# Col 1..6: The 6 CINeMA Domains
# Col 7: Overall Confidence Badge
domain_cols <- c(
  "Within_Study_Bias" = "Within-Study\nBias",
  "Reporting_Bias"    = "Reporting\nBias",
  "Indirectness"      = "Indirectness",
  "Imprecision"       = "Imprecision\n(MCID 0.80-1.25)",
  "Heterogeneity"     = "Heterogeneity\n(\u03c4\u00b2 = 0.009)",
  "Incoherence"       = "Incoherence\n(Node-split)"
)

# Order comparisons from 1 to 15 (reversed for ggplot top-to-bottom layout)
df_cin_ordered <- df_cin[order(df_cin$Comparison_ID, decreasing = TRUE), ]
y_levels <- sprintf("%02d. %s  |  HR: %s", 
                    df_cin_ordered$Comparison_ID, 
                    df_cin_ordered$Comparison, 
                    df_cin_ordered$CI_95)

plot_items <- list()
for (i in 1:nrow(df_cin_ordered)) {
  row <- df_cin_ordered[i, ]
  y_lbl <- y_levels[i]
  
  for (d_idx in seq_along(domain_cols)) {
    d_var <- names(domain_cols)[d_idx]
    d_lbl <- domain_cols[d_idx]
    val <- row[[d_var]]
    
    sym <- switch(val,
      "No concerns"    = "\u2713", # Checkmark / Plus
      "Some concerns"  = "?",
      "Major concerns" = "!",
      "?"
    )
    
    plot_items[[length(plot_items) + 1]] <- data.frame(
      Y_Label = y_lbl,
      Y_Index = i,
      Domain_Var = d_var,
      Col_Index = d_idx,
      Col_Name = d_lbl,
      Rating = val,
      Symbol = sym,
      Is_Confidence = FALSE,
      stringsAsFactors = FALSE
    )
  }
  
  # Overall Confidence Badge (Column 7)
  plot_items[[length(plot_items) + 1]] <- data.frame(
    Y_Label = y_lbl,
    Y_Index = i,
    Domain_Var = "Confidence",
    Col_Index = 7,
    Col_Name = "Overall\nConfidence",
    Rating = row$Confidence,
    Symbol = row$Confidence,
    Is_Confidence = TRUE,
    stringsAsFactors = FALSE
  )
}

df_plot_cin <- do.call(rbind, plot_items)
df_plot_cin$Y_Label <- factor(df_plot_cin$Y_Label, levels = y_levels)

all_col_names <- c(unname(domain_cols), "Overall\nConfidence")
df_plot_cin$Col_Name <- factor(df_plot_cin$Col_Name, levels = all_col_names)

# 3. Aesthetics & Palettes
colors_domain <- c(
  "No concerns"    = "#38A169", # Emerald Green
  "Some concerns"  = "#ECC94B", # Amber / Yellow
  "Major concerns" = "#E53E3E"  # Crimson Red
)

colors_conf <- c(
  "High"      = "#065F46", # Deep Forest Green
  "Moderate"  = "#D97706", # Warm Ochre / Amber
  "Low"       = "#DC2626", # Crimson Red
  "Very Low"  = "#7F1D1D"  # Dark Ruby
)

# 4. Render Publication CINeMA Traffic Light Matrix (300 DPI)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
output_fig <- "outputs/figures/10_cinema_traffic_light.png"
cat(sprintf(" - Rendering Figure 10 to: %s ...\n", output_fig))

p_cin <- ggplot() +
  # Alternating row background
  geom_tile(
    data = df_plot_cin,
    aes(x = Col_Name, y = Y_Label),
    fill = ifelse(as.numeric(df_plot_cin$Y_Label) %% 2 == 0, "#F8FAFC", "#FFFFFF"),
    width = 1, height = 1
  ) +
  # Divider before Overall Confidence
  geom_vline(xintercept = 6.5, color = "#CBD5E0", linewidth = 1.1, linetype = "solid") +
  # Grid lines
  geom_hline(yintercept = seq(0.5, length(y_levels) + 0.5, 1), color = "#EDF2F7", linewidth = 0.5) +
  geom_vline(xintercept = seq(0.5, 7.5, 1), color = "#EDF2F7", linewidth = 0.5) +
  # 6 Domain Traffic Lights (Points)
  geom_point(
    data = subset(df_plot_cin, !Is_Confidence),
    aes(x = Col_Name, y = Y_Label, color = Rating),
    size = 9.0
  ) +
  scale_color_manual(
    name = "Domain Rating:",
    values = colors_domain,
    breaks = c("No concerns", "Some concerns", "Major concerns")
  ) +
  # Domain Math / Text Symbols
  geom_text(
    data = subset(df_plot_cin, !Is_Confidence & Rating == "No concerns"),
    aes(x = Col_Name, y = Y_Label),
    label = "+", color = "#FFFFFF", fontface = "bold", size = 5.0
  ) +
  geom_text(
    data = subset(df_plot_cin, !Is_Confidence & Rating == "Some concerns"),
    aes(x = Col_Name, y = Y_Label),
    label = "?", color = "#1A202C", fontface = "bold", size = 4.2
  ) +
  geom_text(
    data = subset(df_plot_cin, !Is_Confidence & Rating == "Major concerns"),
    aes(x = Col_Name, y = Y_Label),
    label = "-", color = "#FFFFFF", fontface = "bold", size = 5.0
  ) +
  # Overall Confidence Badges (Col 7)
  geom_tile(
    data = subset(df_plot_cin, Is_Confidence),
    aes(x = Col_Name, y = Y_Label, fill = Rating),
    width = 0.88, height = 0.65, color = "#FFFFFF", linewidth = 0.5
  ) +
  scale_fill_manual(
    name = "Overall Confidence:",
    values = colors_conf,
    breaks = c("High", "Moderate", "Low", "Very Low")
  ) +
  geom_text(
    data = subset(df_plot_cin, Is_Confidence),
    aes(x = Col_Name, y = Y_Label, label = Rating),
    color = "#FFFFFF", fontface = "bold", size = 3.6
  ) +
  scale_x_discrete(position = "top") +
  labs(
    title = "CINeMA: Confidence in Network Meta-Analysis Methodological Matrix",
    subtitle = "Systematic Quality & Certainty Assessment Across 15 Pairwise Head-to-Head and Indirect Comparisons",
    x = NULL,
    y = NULL,
    caption = paste0(
      "Framework: CINeMA (Salanti et al., Nikolakopoulou et al. PLOS Med 2020).\n",
      "Minimal Clinically Important Difference (MCID) equivalence range: Hazard Ratio 0.80 to 1.25.\n",
      "Symbols: (+) No concerns; (?) Some concerns; (-) Major concerns."
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15.5, color = "#1A365D", margin = margin(b = 6)),
    plot.subtitle = element_text(color = "#4A5568", size = 11, margin = margin(b = 16)),
    plot.caption = element_text(color = "#718096", size = 9, hjust = 0, lineheight = 1.3, margin = margin(t = 12)),
    axis.text.x.top = element_text(face = "bold", size = 10.5, color = "#2D3748", lineheight = 1.15),
    axis.text.y = element_text(face = "bold", size = 10, color = "#1A202C", family = "sans"),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold", size = 10.5, color = "#1A365D"),
    legend.text = element_text(size = 10, color = "#2D3748"),
    legend.background = element_rect(fill = "#F7FAFC", color = "#E2E8F0", linewidth = 0.5),
    legend.margin = margin(t = 5, b = 5, l = 10, r = 10),
    plot.margin = margin(20, 24, 16, 24)
  )

ggsave(output_fig, plot = p_cin, width = 12.8, height = 9.8, dpi = 300)

cat(sprintf(" [SUCCESS] Figure 10 (CINeMA Traffic Light Matrix) rendered cleanly: %s\n\n", output_fig))
