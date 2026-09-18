# ==============================================================================
# Design Script: scripts/designs/fig06_funnel_plot.R
# Visual Target: Figure 6 - Comparison-Adjusted Funnel Plot (Small-Study Effects)
# Output File:   outputs/figures/06_funnel_plot.png (300 DPI Publication Figure)
# Framework:     netmeta (Chaimani & Salanti Methodology)
# ==============================================================================

suppressPackageStartupMessages({
  library(netmeta)
})

cat("\n======================================================================\n")
cat(" [DESIGN 6/7] FIGURE 6: COMPARISON-ADJUSTED FUNNEL PLOT\n")
cat("======================================================================\n")

# 1. Load Clinical Trial Contrast Data
data_path <- "data/nsclc_trial_contrasts.csv"
if (!file.exists(data_path)) {
  stop(sprintf("Data file not found at: %s. Please run scripts/02_generate_data.R first.", data_path))
}
dat <- read.csv(data_path, stringsAsFactors = FALSE)

# 2. Fit Model
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

# 3. Hierarchy Ordering for Comparison Centering
rk <- netrank(nma, small.values = "good")
trt_order <- names(sort(rk$ranking.random, decreasing = TRUE))

# 4. Render Publication Funnel Plot (300 DPI)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
output_fig <- "outputs/figures/06_funnel_plot.png"
cat(sprintf(" - Rendering Figure 6 to: %s ...\n", output_fig))

png(output_fig, width = 2800, height = 2400, res = 300)
par(mar = c(4.5, 4.5, 3.5, 2))

funnel(
  nma,
  order = trt_order,
  pooled = "random",
  pch = 19,
  col = "#2B6CB0",
  cex = 1.3,
  linreg = TRUE,
  main = "Comparison-Adjusted Funnel Plot (Evaluation of Small-Study Effects)",
  xlab = "Log Hazard Ratio centered by comparison-specific effect",
  ylab = "Standard Error of Log Hazard Ratio"
)

dev.off()

cat(sprintf(" [SUCCESS] Figure 6 rendered cleanly: %s\n\n", output_fig))
