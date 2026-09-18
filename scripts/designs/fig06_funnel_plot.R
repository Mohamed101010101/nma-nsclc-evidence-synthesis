# ==============================================================================
# Design Script: scripts/designs/fig06_funnel_plot.R
# Visual Target: Figure 6 - Comparison-Adjusted Funnel Plot (Small-Study Effects)
# Output File:   outputs/figures/06_funnel_plot.png (300 DPI Publication Figure)
# Framework:     netmeta (Uses Cached NMA Model & Hierarchy)
# ==============================================================================

suppressPackageStartupMessages({
  library(netmeta)
})

cat("\n======================================================================\n")
cat(" [DESIGN 6/7] FIGURE 6: COMPARISON-ADJUSTED FUNNEL PLOT\n")
cat("======================================================================\n")

# 1. Load Cached Model & Rankings (Auto-fit if missing or data changed)
model_path   <- "outputs/models/nma_model.rds"
ranking_path <- "outputs/models/nma_rankings.rds"
data_path    <- "data/nsclc_trial_contrasts.csv"

needs_refit <- !file.exists(model_path) || !file.exists(ranking_path) ||
               (file.exists(data_path) && file.mtime(data_path) > file.mtime(model_path))

if (needs_refit) {
  cat(" - Dataset updated or cached model missing. Re-fitting model via 01_fit_nma_model.R ...\n")
  refit_env <- new.env(parent = globalenv())
  refit_env$force_refit <- TRUE
  source("scripts/analyses/01_fit_nma_model.R", local = refit_env)
}

nma <- readRDS(model_path)
rk <- readRDS(ranking_path)
trt_order <- names(sort(rk$ranking.random, decreasing = TRUE))
cat(" - Loaded cached model & rankings in < 0.01 seconds.\n")

# 2. Render Publication Funnel Plot (300 DPI)
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
