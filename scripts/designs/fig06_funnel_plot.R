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
cat(" [DESIGN 06/12] FIGURE 06: COMPARISON-ADJUSTED FUNNEL PLOT\n")
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

load_start_time <- Sys.time()
nma <- readRDS(model_path)
rk <- readRDS(ranking_path)
trt_order <- names(sort(rk$ranking.random, decreasing = TRUE))
load_end_time <- Sys.time()
load_time_taken <- round(as.numeric(difftime(load_end_time, load_start_time, units="secs")), 3)
cat(sprintf(" - Loaded cached model & rankings in %.3f seconds.\n", load_time_taken))

trt_labels_map <- c(
  "IO_Chemo"  = "IO + Chemo",
  "TKI_Chemo" = "TKI + Chemo",
  "Dual_IO"   = "Dual IO",
  "IO_Mono"   = "IO Monotherapy",
  "TKI"       = "TKI Monotherapy",
  "Chemo"     = "Chemotherapy"
)

# 2. Render Publication Funnel Plot (300 DPI)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
output_fig <- "outputs/figures/06_funnel_plot.png"
cat(sprintf(" - Rendering Figure 6 to: %s ...\n", output_fig))

png(output_fig, width = 2800, height = 2400, res = 300)
par(mar = c(4.5, 4.5, 3.5, 2))

fn_res <- funnel(
  nma,
  order = trt_order,
  pooled = "random",
  pch = 19,
  col = "#2B6CB0",
  cex = 1.3,
  linreg = TRUE,
  legend = FALSE,
  main = "Comparison-Adjusted Funnel Plot (Evaluation of Small-Study Effects)",
  xlab = "Log Hazard Ratio centered by comparison-specific effect",
  ylab = "Standard Error of Log Hazard Ratio"
)

# Clean comparison labels for publication legend
raw_comps <- unique(fn_res$comparison)
clean_comps <- sapply(raw_comps, function(comp) {
  parts <- strsplit(comp, ":")[[1]]
  p1 <- ifelse(parts[1] %in% names(trt_labels_map), trt_labels_map[parts[1]], parts[1])
  p2 <- ifelse(parts[2] %in% names(trt_labels_map), trt_labels_map[parts[2]], parts[2])
  paste(p1, "vs", p2)
})

legend(
  "topright",
  legend = clean_comps,
  pch = 19,
  col = "#2B6CB0",
  bty = "o",
  box.col = "#CBD5E0",
  bg = "#FFFFFFEE",
  cex = 0.85
)

dev.off()

cat(sprintf(" [SUCCESS] Figure 6 rendered cleanly: %s\n\n", output_fig))
