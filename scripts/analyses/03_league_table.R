# ==============================================================================
# Script: scripts/analyses/03_league_table.R
# Purpose: Pairwise Dual-Model League Table Generation (CSV Matrix)
# Inputs:  outputs/models/nma_model.rds & outputs/models/nma_rankings.rds
# Output:  outputs/tables/league_table_random_common.csv
# Package: netmeta (Rücker & Schwarzer)
# ==============================================================================

suppressPackageStartupMessages({
  library(netmeta)
})

cat("\n======================================================================\n")
cat(" [ANALYSIS 3/5] LEAGUE TABLE GENERATION (CSV MATRIX)\n")
cat("======================================================================\n")

# 1. Load Cached Model & Rankings (Auto-fit if missing)
model_path <- "outputs/models/nma_model.rds"
ranking_path <- "outputs/models/nma_rankings.rds"

if (!file.exists(model_path) || !file.exists(ranking_path)) {
  cat(" - Cached model not detected. Running scripts/analyses/01_fit_nma_model.R ...\n")
  source("scripts/analyses/01_fit_nma_model.R", local = new.env())
}

nma <- readRDS(model_path)
rk <- readRDS(ranking_path)
trt_order <- names(sort(rk$ranking.random, decreasing = TRUE))
cat(" - Loaded cached model & rankings in < 0.01 seconds.\n")

# 2. Construct Dual-Model League Table
# Lower triangle: Random-effects model HR [95% CI]
# Upper triangle: Common-effects model HR [95% CI]
lg <- netleague(nma, digits = 2, seq = trt_order)

# 3. Export Raw CSV Matrix
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
output_csv <- "outputs/tables/league_table_random_common.csv"
write.csv(lg$random, output_csv)
cat(sprintf(" - Exported raw league matrix CSV to: %s\n\n", output_csv))
print(lg$random)
cat("\n")
