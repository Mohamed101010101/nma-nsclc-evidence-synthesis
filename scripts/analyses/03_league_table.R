# ==============================================================================
# Script: scripts/analyses/03_league_table.R
# Purpose: Pairwise Dual-Model League Table Generation (CSV Table Output)
# Output:  outputs/tables/league_table_random_common.csv
# Package: netmeta (Rücker & Schwarzer)
# ==============================================================================

suppressPackageStartupMessages({
  library(netmeta)
})

cat("\n======================================================================\n")
cat(" [ANALYSIS 3/5] LEAGUE TABLE GENERATION (CSV MATRIX)\n")
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

# 3. Determine Ordering by Hierarchy (P-scores)
rk <- netrank(nma, small.values = "good")
pscores_rand <- rk$ranking.random
trt_order <- names(sort(pscores_rand, decreasing = TRUE))

# 4. Construct Dual-Model League Table
# Lower triangle: Random-effects model HR [95% CI]
# Upper triangle: Common-effects model HR [95% CI]
lg <- netleague(nma, digits = 2, seq = trt_order)

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
output_csv <- "outputs/tables/league_table_random_common.csv"
write.csv(lg$random, output_csv)
cat(sprintf(" - Exported raw league matrix CSV to: %s\n\n", output_csv))
print(lg$random)
cat("\n")
