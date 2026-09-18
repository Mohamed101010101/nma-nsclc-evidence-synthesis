# ==============================================================================
# Script: 04_forest_plot.R
# Purpose: Forest Plot of Active Treatment Regimens vs Reference Control (Chemo)
# Output:  outputs/figures/02_forest_plot_random.png (300 DPI Publication Figure)
# Package: netmeta (Frequentist Random-Effects Model)
# ==============================================================================

suppressPackageStartupMessages({
  library(netmeta)
})

cat("\n======================================================================\n")
cat(" [ANALYSIS 2/7] FOREST PLOT: ACTIVE REGIMENS VS CHEMOTHERAPY\n")
cat("======================================================================\n")

# 1. Load Clinical Trial Contrast Data
data_path <- "data/nsclc_trial_contrasts.csv"
if (!file.exists(data_path)) {
  stop(sprintf("Data file not found at: %s. Please run scripts/02_generate_data.R first.", data_path))
}
dat <- read.csv(data_path, stringsAsFactors = FALSE)
cat(sprintf(" - Loaded contrast dataset: %d comparisons across %d trials\n",
            nrow(dat), length(unique(dat$studlab))))

# 2. Fit Frequentist Graph-Theoretical Model
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

# 3. Compute Treatment Rankings for Hierarchical Forest Plot Ordering
rk <- netrank(nma, small.values = "good")

# 4. Render Publication Forest Plot (300 DPI)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
output_fig <- "outputs/figures/02_forest_plot_random.png"
cat(sprintf(" - Rendering Figure 2 to: %s ...\n", output_fig))

png(output_fig, width = 3200, height = 1800, res = 300)

forest(
  nma,
  reference.group = "Chemo",
  pooled = "random",
  sortvar = -rk$ranking.random,
  smlab = "Hazard Ratio (95% CI)\nvs Chemotherapy",
  label.left = "Favors Active Regimen",
  label.right = "Favors Chemotherapy",
  drop.reference.group = TRUE,
  digits = 2,
  col.square = "#2B6CB0",
  col.diamond = "#C53030",
  col.inside = "#1A202C",
  header.line = TRUE,
  leftcols = c("studlab"),
  leftlabs = c("Treatment Regimen"),
  rightcols = c("effect", "ci"),
  rightlabs = c("HR", "95% CI")
)

dev.off()

cat(sprintf(" [SUCCESS] Forest Plot saved cleanly: %s\n\n", output_fig))
