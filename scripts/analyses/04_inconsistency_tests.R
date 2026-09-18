# ==============================================================================
# Script: scripts/analyses/04_inconsistency_tests.R
# Purpose: Global Q Variance Decomposition & Local Inconsistency Diagnostics
# Output:  outputs/tables/inconsistency_statistics.csv
# Package: netmeta (Design-by-Treatment Interaction Model)
# ==============================================================================

suppressPackageStartupMessages({
  library(netmeta)
})

cat("\n======================================================================\n")
cat(" [ANALYSIS 4/4] INCONSISTENCY EVALUATION (GLOBAL Q & NODE-SPLITTING)\n")
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

# 3. Global Inconsistency: Decomposition of Cochran's Q
# Q_total = Q_within (heterogeneity) + Q_between (inconsistency)
df_inconsistency <- data.frame(
  Source = c("Total Variation (Q)", "Within-Designs Heterogeneity (Q_het)", "Between-Designs Inconsistency (Q_inc)"),
  Q_Statistic = round(c(nma$Q, nma$Q.heterogeneity, nma$Q.inconsistency), 2),
  Degrees_of_Freedom = c(nma$df.Q, nma$df.Q.heterogeneity, nma$df.Q.inconsistency),
  P_Value = sprintf("%.4f", c(nma$pval.Q, nma$pval.Q.heterogeneity, nma$pval.Q.inconsistency)),
  Interpretation = c(
    ifelse(nma$pval.Q > 0.05, "No significant total excess variance", "Significant total variation"),
    ifelse(nma$pval.Q.heterogeneity > 0.05, "Homogeneity within trial designs", "Heterogeneity within designs"),
    ifelse(nma$pval.Q.inconsistency > 0.05, "Full Transitivity/Consistency upheld", "Evidence of Inconsistency")
  ),
  stringsAsFactors = FALSE
)

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
output_tbl <- "outputs/tables/inconsistency_statistics.csv"
write.csv(df_inconsistency, output_tbl, row.names = FALSE)
cat(sprintf(" - Exported global inconsistency table to: %s\n\n", output_tbl))
print(df_inconsistency)
cat("\n")
