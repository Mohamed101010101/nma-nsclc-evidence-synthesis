# ==============================================================================
# Script: scripts/analyses/01_fit_nma_model.R
# Purpose: Core Frequentist Graph-Theoretical Network Meta-Analysis Model Fitting
# Package: netmeta (Rücker 2012 Electrical Network Analogy)
# Output:  Model summary, heterogeneity diagnostics (tau2, I2)
# ==============================================================================

suppressPackageStartupMessages({
  library(netmeta)
})

cat("\n======================================================================\n")
cat(" [ANALYSIS 1/4] MODEL ESTIMATION & HETEROGENEITY ASSESSMENT\n")
cat("======================================================================\n")

# 1. Load Clinical Trial Contrast Data
data_path <- "data/nsclc_trial_contrasts.csv"
if (!file.exists(data_path)) {
  stop(sprintf("Data file not found at: %s. Please run scripts/02_generate_data.R first.", data_path))
}
dat <- read.csv(data_path, stringsAsFactors = FALSE)
cat(sprintf(" - Loaded contrast dataset: %d comparisons across %d trials\n",
            nrow(dat), length(unique(dat$studlab))))
cat(sprintf(" - Total patients evaluated: %s\n", 
            format(sum(dat$n_total[!duplicated(dat$studlab)]), big.mark = ",")))

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

cat("\n [MODEL SUMMARY]\n")
cat(sprintf(" - Treatments (n): %d\n", nma$n))
cat(sprintf(" - Pairwise comparisons (m): %d\n", nma$m))
cat(sprintf(" - Study designs (d): %d\n", nma$d))
cat(sprintf(" - Between-study heterogeneity (tau^2): %.4f (tau = %.4f)\n", nma$tau2, nma$tau))
cat(sprintf(" - Inconsistency/Heterogeneity index (I^2): %.1f%%\n", nma$I2 * 100))
cat(sprintf(" - Total Cochran's Q: %.2f (df = %d, p = %.4f)\n", nma$Q, nma$df.Q, nma$pval.Q))

cat("\n [SUCCESS] Network Meta-Analysis model fitted successfully.\n\n")
