# ==============================================================================
# Design Script: scripts/designs/fig04_netsplit_inconsistency.R
# Visual Target: Figure 4 - Node-Splitting Local Inconsistency Forest Plot
# Output File:   outputs/figures/04_netsplit_inconsistency.png (300 DPI Publication Figure)
# Framework:     netmeta (Node-Splitting / Separate Direct vs Indirect Evidence)
# ==============================================================================

suppressPackageStartupMessages({
  library(netmeta)
})

cat("\n======================================================================\n")
cat(" [DESIGN 4/7] FIGURE 4: NODE-SPLITTING LOCAL INCONSISTENCY FOREST PLOT\n")
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

# 3. Calculate Node-Splitting Models
cat(" - Computing node-splitting models across all closed evidence loops ...\n")
ns <- netsplit(nma)

# 4. Render Publication Node-Splitting Forest Plot (300 DPI)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
output_fig <- "outputs/figures/04_netsplit_inconsistency.png"
cat(sprintf(" - Rendering Figure 4 to: %s ...\n", output_fig))

# Full dimensions (3400x4600px) ensure zero truncation of loops, tests, and axis
png(output_fig, width = 3400, height = 4600, res = 300)

forest(
  ns,
  pooled = "random",
  fontsize = 9,
  spacing = 1.05,
  digits = 2,
  smlab = "Hazard Ratio (95% CI)\nDirect vs Indirect vs Network"
)

dev.off()

cat(sprintf(" [SUCCESS] Figure 4 rendered cleanly: %s\n\n", output_fig))
