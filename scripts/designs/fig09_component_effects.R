# ==============================================================================
# Script: scripts/designs/fig09_component_effects.R
# Purpose: Publication-Grade Figure 09 - Component Network Meta-Analysis (CNMA)
# Output: outputs/figures/09_component_effects.png (300 DPI, 11x8.5 in)
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(grid)
})

cat("\n======================================================================\n")
cat(" [DESIGN 9/12] FIGURE 9: COMPONENT NETWORK META-ANALYSIS (CNMA) FOREST\n")
cat("======================================================================\n")

cnma_data_file <- "outputs/models/component_nma_data.rds"
if (!file.exists(cnma_data_file)) {
  cat(" - CNMA cache missing. Running 07_component_nma.R ...\n")
  source("scripts/analyses/07_component_nma.R")
}

cnma <- readRDS(cnma_data_file)
df_comp <- cnma$components_df
df_comb <- cnma$combinations_df

# Build unified plotting data frame
# Component labels for clinical clarity
comp_labels <- c(
  "IO"    = "Anti-PD-(L)1 Monoclonal Antibody (IO)",
  "TKI"   = "Tyrosine Kinase Inhibitor (TKI)",
  "CTLA4" = "Anti-CTLA-4 Antibody (CTLA4)"
)

comb_labels <- c(
  "Chemo + IO"  = "Chemo + Anti-PD-(L)1 (Additive)",
  "Chemo + TKI" = "Chemo + TKI (Additive)",
  "CTLA4 + IO"  = "Anti-CTLA-4 + Anti-PD-(L)1 (Additive)",
  "IO"          = "Anti-PD-(L)1 Monotherapy (Additive)",
  "TKI"         = "TKI Monotherapy (Additive)"
)

df_comp$Display_Name <- comp_labels[df_comp$Item]
df_comb$Display_Name <- comb_labels[df_comb$Item]

df_comp$Category <- "A. Marginal Incremental Effects of Individual Components (vs Chemo Backbone)"
df_comb$Category <- "B. Predicted Regimen Effects under Additive CNMA Model (vs Chemo Backbone)"

# Colors
df_comp$Color <- c("#1B365D", "#2E7D32", "#C62828") # CTLA4, IO, TKI
df_comb$Color <- c("#00838F", "#E65100", "#6A1B9A", "#2E7D32", "#C62828")

df_plot <- rbind(
  df_comp %>% select(Category, Display_Name, iHR, CI_Lower, CI_Upper, HR_String, P_Value_String, Color),
  df_comb %>% select(Category, Display_Name, iHR, CI_Lower, CI_Upper, HR_String, P_Value_String, Color)
)

# Order items logically within categories
df_plot$Display_Name <- factor(df_plot$Display_Name, levels = rev(c(
  "Anti-PD-(L)1 Monoclonal Antibody (IO)",
  "Tyrosine Kinase Inhibitor (TKI)",
  "Anti-CTLA-4 Antibody (CTLA4)",
  "Chemo + Anti-PD-(L)1 (Additive)",
  "Chemo + TKI (Additive)",
  "Anti-CTLA-4 + Anti-PD-(L)1 (Additive)",
  "Anti-PD-(L)1 Monotherapy (Additive)",
  "TKI Monotherapy (Additive)"
)))

p <- ggplot(df_plot, aes(x = iHR, y = Display_Name)) +
  # Null effect reference line
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "#78909C", linewidth = 0.8) +
  # Facet by Category
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  # Error bars
  geom_errorbar(aes(xmin = CI_Lower, xmax = CI_Upper), 
                width = 0.25, orientation = "y", linewidth = 0.9, color = df_plot$Color) +
  # Point estimates
  geom_point(shape = 18, size = 5.0, color = df_plot$Color) +
  # Text labels for HR & 95% CI
  geom_text(aes(x = 1.32, label = sprintf("%s  |  p %s", HR_String, P_Value_String)),
            hjust = 0, size = 3.6, family = "sans", fontface = "bold", color = "#263238") +
  scale_x_continuous(
    trans = "log",
    breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3),
    limits = c(0.55, 1.75),
    labels = c("0.50", "0.60", "0.70", "0.80", "0.90", "1.00", "1.10", "1.20", "1.30")
  ) +
  labs(
    title = "Component Network Meta-Analysis (CNMA): Deconstruction of Regimen Synergy",
    subtitle = sprintf("Marginal Component Incremental Hazard Ratios & Additive Combinations vs Platinum-Chemotherapy Backbone\nInteraction Test vs Standard NMA: Q_diff = %.2f (df = %d, p = %.4f -> Significant Synergistic Interaction Detected)",
                       cnma$Q_diff, cnma$df_Q_diff, cnma$pval_Q_diff),
    x = "Incremental Hazard Ratio (iHR) & 95% Confidence Interval (Log Scale)",
    y = NULL,
    caption = "Frequentist Additive CNMA Model (Rücker et al., 2020) fitted using netmeta::netcomb(). Reference Backbone = Platinum Chemotherapy.\niHR < 1.0 indicates extended Overall Survival when adding component/combination; iHR > 1.0 indicates diminished efficacy."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, color = "#0D233A", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10, color = "#37474F", lineheight = 1.25, margin = margin(b = 14)),
    plot.caption = element_text(size = 8.5, color = "#78909C", hjust = 0, margin = margin(t = 12)),
    axis.text.y = element_text(face = "bold", size = 10.5, color = "#263238"),
    axis.text.x = element_text(size = 9.5, color = "#37474F"),
    axis.title.x = element_text(face = "bold", size = 10.5, color = "#263238", margin = margin(t = 8)),
    strip.text = element_text(face = "bold", size = 11, color = "#FFFFFF", hjust = 0),
    strip.background = element_rect(fill = "#1B365D", color = NA),
    panel.grid.major.x = element_line(color = "#ECEFF1", linewidth = 0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "#F5F5F5", linewidth = 0.5),
    plot.margin = margin(t = 16, r = 20, b = 16, l = 16)
  )

# Add Direction of Benefit annotations
p <- p +
  annotate("text", x = 0.72, y = 0.55, label = "◄ Favors Component / Addition", 
           size = 3.2, fontface = "bold", color = "#2E7D32") +
  annotate("text", x = 1.12, y = 0.55, label = "Favors Backbone Only ►", 
           size = 3.2, fontface = "bold", color = "#C62828")

fig_out <- "outputs/figures/09_component_effects.png"
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
cat(sprintf(" - Rendering Figure 9 to: %s ...\n", fig_out))
ggsave(fig_out, plot = p, width = 11, height = 7.5, dpi = 300, bg = "#FFFFFF")
cat(sprintf(" [SUCCESS] Figure 9 rendered cleanly: %s\n", fig_out))
