# ==============================================================================
# Script: 07_inconsistency_netsplit.R
# Purpose: Global Inconsistency Decomposition & Local Node-Splitting Analysis
# Outputs: outputs/tables/inconsistency_statistics.csv
#          outputs/figures/04_netsplit_inconsistency.png (300 DPI Publication Figure)
# Package: netmeta (Frequentist Design-by-Treatment Interaction Model)
# ==============================================================================

suppressPackageStartupMessages({
  library(netmeta)
})

cat("\n======================================================================\n")
cat(" [ANALYSIS 5/7] INCONSISTENCY EVALUATION: GLOBAL Q & NODE-SPLITTING\n")
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
cat(sprintf(" - Exported global inconsistency table to: %s\n", output_tbl))
print(df_inconsistency)

# 4. Local Inconsistency via Node-Splitting
# Separates direct head-to-head evidence from indirect evidence for each closed loop
cat("\n - Calculating node-splitting models for all closed loops ...\n")
ns <- netsplit(nma)

# 5. Render Publication Node-Splitting Forest Plot (300 DPI)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
output_fig <- "outputs/figures/04_netsplit_inconsistency.png"
cat(sprintf(" - Rendering Figure 4 to: %s ...\n", output_fig))

# Dimensions (3400x4600px) ensure all comparisons and axis values are fully visible
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

cat(sprintf(" [SUCCESS] Node-splitting analysis & Figure 4 saved cleanly: %s\n\n", output_fig))
