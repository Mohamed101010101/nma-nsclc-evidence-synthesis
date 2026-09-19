# ==============================================================================
# Script: scripts/designs/fig11_benefit_risk_tradeoff.R
# Purpose: Publication-Grade Figure 11 - Bi-dimensional Benefit-Risk Trade-Off Matrix
# Output: outputs/figures/11_benefit_risk_tradeoff.png (300 DPI, 11x8.5 in)
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(grid)
})

cat("\n======================================================================\n")
cat(" [DESIGN 11/12] FIGURE 11: BI-DIMENSIONAL BENEFIT-RISK TRADE-OFF MATRIX\n")
cat("======================================================================\n")

br_data_file <- "outputs/models/benefit_risk_data.rds"
if (!file.exists(br_data_file)) {
  cat(" - Benefit-risk data cache missing. Running 09_benefit_risk_tradeoff.R ...\n")
  source("scripts/analyses/09_benefit_risk_tradeoff.R")
}

br_data <- readRDS(br_data_file)
df_br   <- br_data$benefit_risk_df

# Regimen labels and colors
trt_labels <- c(
  "IO_Chemo"  = "IO + Chemo",
  "TKI_Chemo" = "TKI + Chemo",
  "Dual_IO"   = "Dual IO",
  "IO_Mono"   = "IO Monotherapy",
  "TKI"       = "TKI Monotherapy",
  "Chemo"     = "Chemo Alone (Ref)"
)

palette_regimens <- c(
  "IO_Chemo"  = "#1B365D",  # Navy
  "TKI_Chemo" = "#E65100",  # Orange
  "Dual_IO"   = "#6A1B9A",  # Purple
  "IO_Mono"   = "#2E7D32",  # Forest Green
  "TKI"       = "#C62828",  # Red
  "Chemo"     = "#546E7A"   # Grey
)

df_br$Label <- trt_labels[df_br$Treatment]
df_br$Color <- palette_regimens[df_br$Treatment]

# Invert Y-axis conceptually: we want superior survival (lower HR) at the top!
# In ggplot, scale_y_reverse() achieves this naturally!

p <- ggplot(df_br, aes(x = OR_Tox, y = HR_OS)) +
  # Quadrant Shading
  # Quadrant II (Top-Left): High Efficacy, Low Toxicity
  annotate("rect", xmin = 0.20, xmax = 1.00, ymin = 0.60, ymax = 0.82,
           fill = "#E8F5E9", alpha = 0.55) +
  # Quadrant I (Top-Right): High Efficacy, High Toxicity
  annotate("rect", xmin = 1.00, xmax = 2.40, ymin = 0.60, ymax = 0.82,
           fill = "#EDE7F6", alpha = 0.55) +
  # Quadrant III (Bottom-Left): Modest Efficacy, Low Toxicity
  annotate("rect", xmin = 0.20, xmax = 1.00, ymin = 0.82, ymax = 1.08,
           fill = "#FFF8E1", alpha = 0.55) +
  # Quadrant IV (Bottom-Right): Modest Efficacy, High Toxicity
  annotate("rect", xmin = 1.00, xmax = 2.40, ymin = 0.82, ymax = 1.08,
           fill = "#FFEBEE", alpha = 0.55) +
  
  # Quadrant Header Labels
  annotate("text", x = 0.23, y = 0.615, hjust = 0, vjust = 1,
           label = "QUADRANT II: OPTIMAL WINDOW\nSuperior Survival + Lower Toxicity",
           size = 3.3, fontface = "bold", color = "#2E7D32", lineheight = 1.1) +
  annotate("text", x = 2.37, y = 0.615, hjust = 1, vjust = 1,
           label = "QUADRANT I: INTENSIVE COMBINATION\nSuperior Survival + Increased Toxicity",
           size = 3.3, fontface = "bold", color = "#4527A0", lineheight = 1.1) +
  annotate("text", x = 0.23, y = 1.065, hjust = 0, vjust = 0,
           label = "QUADRANT III: TOLERABLE COMPROMISE\nModest Survival + Low Severe Toxicity",
           size = 3.3, fontface = "bold", color = "#F57F17", lineheight = 1.1) +
  annotate("text", x = 2.37, y = 1.065, hjust = 1, vjust = 0,
           label = "QUADRANT IV: UNFAVORABLE BACKBONE\nStandard Survival + Standard/High Toxicity",
           size = 3.3, fontface = "bold", color = "#C62828", lineheight = 1.1) +

  # Reference dashed lines for Chemotherapy Backbone
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "#37474F", linewidth = 0.8) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "#37474F", linewidth = 0.8) +
  
  # Bidirectional 95% Confidence Intervals
  # Horizontal: Toxicity OR 95% CI
  geom_errorbar(aes(xmin = OR_Tox_Lower, xmax = OR_Tox_Upper), 
                orientation = "y", width = 0.015, linewidth = 0.8, color = df_br$Color, alpha = 0.75) +
  # Vertical: Overall Survival HR 95% CI
  geom_errorbar(aes(ymin = HR_OS_Lower, ymax = HR_OS_Upper), 
                orientation = "x", width = 0.04, linewidth = 0.8, color = df_br$Color, alpha = 0.75) +
  
  # Treatment Point Estimates (Size scaled by Net Clinical Benefit)
  geom_point(aes(size = Net_Benefit_Score), color = "#263238", fill = df_br$Color, 
             shape = 21, stroke = 1.5, alpha = 0.95) +
  scale_size_continuous(range = c(5.5, 9.5), guide = "none") +
  
  # Repositioned Text Annotations to prevent overlap
  geom_label(aes(label = sprintf("%s\nHR: %.2f | OR: %.2f", Label, HR_OS, OR_Tox)),
             vjust = c(-0.5, 1.4, -0.6, 1.3, -0.5, 1.4)[match(df_br$Treatment, c("IO_Chemo", "TKI_Chemo", "Dual_IO", "IO_Mono", "TKI", "Chemo"))],
             hjust = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)[match(df_br$Treatment, c("IO_Chemo", "TKI_Chemo", "Dual_IO", "IO_Mono", "TKI", "Chemo"))],
             size = 3.2, fontface = "bold", fill = "#FFFFFF", color = "#263238",
             label.padding = unit(0.2, "lines"), label.r = unit(0.15, "lines")) +
  
  # Axis Scaling
  scale_x_continuous(
    trans = "log",
    breaks = c(0.25, 0.35, 0.50, 0.70, 1.00, 1.40, 1.80, 2.30),
    limits = c(0.20, 2.45),
    labels = c("0.25", "0.35", "0.50", "0.70", "1.00", "1.40", "1.80", "2.30")
  ) +
  scale_y_reverse(
    trans = "log",
    breaks = c(0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05),
    limits = c(1.08, 0.60),
    labels = c("0.60", "0.65", "0.70", "0.75", "0.80", "0.85", "0.90", "0.95", "1.00", "1.05")
  ) +
  labs(
    title = "Bi-dimensional Benefit-Risk Trade-Off Matrix: Survival Efficacy vs Severe Toxicity",
    subtitle = "Simultaneous Dual Network Meta-Analysis Mapping Overall Survival Hazard Ratio against Grade 3-5 Adverse Event Odds Ratio\nReference Arm: Platinum Chemotherapy Backbone (HR = 1.00, OR = 1.00; Point size proportional to Net Clinical Benefit Score)",
    x = "Severe Toxicity (Grade 3-5 AEs): Odds Ratio vs Chemotherapy (Log Scale, < 1.0 = Safer)",
    y = "Overall Survival Efficacy: Hazard Ratio vs Chemotherapy (Inverted Log Scale, Top = Superior Survival)",
    caption = "Derived from dual frequentist random-effects network meta-analyses across 24 randomized controlled trials (N = 14,357 patients).\nHorizontal and vertical error bars denote 95% confidence intervals. Quadrant thresholds define clinical decision domains."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15, color = "#0D233A", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10.5, color = "#37474F", lineheight = 1.25, margin = margin(b = 14)),
    plot.caption = element_text(size = 8.5, color = "#78909C", hjust = 0, margin = margin(t = 12)),
    axis.text = element_text(face = "bold", size = 9.5, color = "#263238"),
    axis.title = element_text(face = "bold", size = 10.5, color = "#263238"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.major = element_line(color = "#CFD8DC", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 16, r = 20, b = 16, l = 16)
  )

fig_out <- "outputs/figures/11_benefit_risk_tradeoff.png"
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
cat(sprintf(" - Rendering Figure 11 to: %s ...\n", fig_out))
ggsave(fig_out, plot = p, width = 11.5, height = 8.0, dpi = 300, bg = "#FFFFFF")
cat(sprintf(" [SUCCESS] Figure 11 rendered cleanly: %s\n", fig_out))
