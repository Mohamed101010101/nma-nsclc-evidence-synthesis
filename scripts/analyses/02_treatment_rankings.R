# ==============================================================================
# Script: scripts/analyses/02_treatment_rankings.R
# Purpose: Treatment Hierarchy & P-Score Calculations (Frequentist SUCRA)
# Inputs:  outputs/models/nma_model.rds & outputs/models/nma_rankings.rds
# Output:  outputs/tables/treatment_rankings.csv
# Package: netmeta
# ==============================================================================

suppressPackageStartupMessages({
  library(netmeta)
})

cat("\n======================================================================\n")
cat(" [ANALYSIS 2/5] TREATMENT RANKING VIA P-SCORES (SUCRA ANALOGUE)\n")
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
cat(" - Successfully loaded cached model & rankings in < 0.01 seconds.\n")

# 2. Extract P-Scores & Treatment Hierarchy
pscores_rand <- rk$ranking.random
trt_order <- names(sort(pscores_rand, decreasing = TRUE))

trt_labels_map <- c(
  "IO_Chemo"  = "IO + Chemo",
  "TKI_Chemo" = "TKI + Chemo",
  "Dual_IO"   = "Dual IO",
  "IO_Mono"   = "IO Monotherapy",
  "TKI"       = "TKI Monotherapy",
  "Chemo"     = "Chemotherapy"
)

# 3. Construct Comprehensive Treatment Ranking Table
df_rankings <- data.frame(
  Treatment = trt_labels_map[trt_order],
  Code = trt_order,
  Rank = 1:length(trt_order),
  Pscore_Random = round(pscores_rand[trt_order], 4),
  Pscore_Common = round(rk$ranking.common[trt_order], 4),
  HR_vs_Chemo_Random = ifelse(trt_order == "Chemo", "1.00 (Reference)",
                              sprintf("%.2f [%.2f; %.2f]", 
                                      exp(nma$TE.random[trt_order, "Chemo"]),
                                      exp(nma$lower.random[trt_order, "Chemo"]),
                                      exp(nma$upper.random[trt_order, "Chemo"]))),
  Pval_vs_Chemo = ifelse(trt_order == "Chemo", "Reference",
                         ifelse(nma$pval.random[trt_order, "Chemo"] < 0.0001, "< 0.0001",
                                sprintf("%.4f", nma$pval.random[trt_order, "Chemo"]))),
  stringsAsFactors = FALSE
)

# 4. Export Ranking Table to CSV
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
output_tbl <- "outputs/tables/treatment_rankings.csv"
write.csv(df_rankings[, c("Treatment", "Rank", "Pscore_Random", "Pscore_Common", "HR_vs_Chemo_Random", "Pval_vs_Chemo")], 
          output_tbl, row.names = FALSE)

cat(sprintf(" - Successfully exported ranking table to: %s\n\n", output_tbl))
print(df_rankings[, c("Treatment", "Rank", "Pscore_Random", "HR_vs_Chemo_Random", "Pval_vs_Chemo")])
cat("\n")
