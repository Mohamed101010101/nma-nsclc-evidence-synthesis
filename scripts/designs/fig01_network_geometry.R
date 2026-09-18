# ==============================================================================
# Design Script: scripts/designs/fig01_network_geometry.R
# Visual Target: Figure 1 - Evidence Network Geometry (Network Topology)
# Output File:   outputs/figures/01_network_geometry.png (300 DPI Publication Figure)
# Framework:     netmeta (Frequentist Graph-Theoretical Framework)
# ==============================================================================

suppressPackageStartupMessages({
  library(netmeta)
})

cat("\n======================================================================\n")
cat(" [DESIGN 1/7] FIGURE 1: EVIDENCE NETWORK GEOMETRY\n")
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

# 3. Node & Palette Configuration
colors_nodes <- c(
  "Chemo"     = "#718096", # Slate Grey (Standard Reference)
  "IO_Mono"   = "#3182CE", # Classic Blue (Active Monotherapy)
  "IO_Chemo"  = "#2B6CB0", # Deep Blue (Chemo Combo)
  "Dual_IO"   = "#805AD5", # Purple (Dual Checkpoint)
  "TKI"       = "#D69E2E", # Warm Amber (Targeted Mono)
  "TKI_Chemo" = "#DD6B20"  # Rust Orange (Targeted Combo)
)

# Compute cumulative patient sample size per treatment node
pts_size <- sapply(nma$trts, function(t) {
  sum(dat$n_treat1[dat$treat1 == t], dat$n_treat2[dat$treat2 == t], na.rm = TRUE)
})
pts_cex <- 6.5 + (pts_size / max(pts_size)) * 4.5

# 4. Render Publication Network Geometry (300 DPI)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
output_fig <- "outputs/figures/01_network_geometry.png"
cat(sprintf(" - Rendering Figure 1 to: %s ...\n", output_fig))

png(output_fig, width = 3000, height = 2600, res = 300)
par(mar = c(5.2, 2.5, 3.8, 2.5))

netgraph(
  nma,
  points = TRUE,
  cex.points = pts_cex,
  col.points = colors_nodes[nma$trts],
  col = "#718096",
  plastic = FALSE,
  thickness = "number.of.studies",
  lwd.max = 7.5,
  lwd.min = 2,
  cex = 1.3,
  offset = 0.045,
  multiarm = TRUE,
  col.multiarm = "#E2E8F0",
  main = "Evidence Network Geometry: First-Line NSCLC Overall Survival"
)
mtext("Node diameter proportional to sample size | Line thickness proportional to trial count", 
      side = 3, line = 0.5, cex = 1.0, col = "#4A5568")

# Horizontal centered legend at bottom with optimized offset
legend("bottom", 
       legend = c("1 Trial", "2-3 Trials", "5+ Trials"), 
       lwd = c(2, 4.5, 7.5), 
       col = "#718096", 
       horiz = TRUE,
       bty = "o", 
       box.col = "#CBD5E0",
       bg = "#FFFFFFEE",
       title = "Direct Evidence Base (Line Thickness)", 
       cex = 0.95,
       inset = c(0, -0.045),
       xpd = TRUE)
dev.off()

cat(sprintf(" [SUCCESS] Figure 1 rendered cleanly: %s\n\n", output_fig))
