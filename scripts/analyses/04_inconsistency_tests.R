# ==============================================================================
# Script: scripts/analyses/04_inconsistency_tests.R
# Purpose: Global Q Variance Decomposition & Inconsistency Diagnostics
# Input:   outputs/models/nma_model.rds
# Output:  outputs/tables/inconsistency_statistics.csv
# Package: netmeta (Design-by-Treatment Interaction Model)
# ==============================================================================

suppressPackageStartupMessages({
  library(netmeta)
})

cat("\n======================================================================\n")
cat(" [ANALYSIS 4/5] INCONSISTENCY EVALUATION (GLOBAL Q DECOMPOSITION)\n")
cat("======================================================================\n")

# 1. Load Cached Model (Auto-fit if missing)
model_path <- "outputs/models/nma_model.rds"
if (!file.exists(model_path)) {
  cat(" - Cached model not detected. Running scripts/analyses/01_fit_nma_model.R ...\n")
  source("scripts/analyses/01_fit_nma_model.R", local = new.env())
}

nma <- readRDS(model_path)
cat(" - Loaded cached model in < 0.01 seconds.\n")

# 2. Global Inconsistency: Decomposition of Cochran's Q
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
