# ==============================================================================
# Script: scripts/analyses/02_treatment_rankings.R
# Purpose: Treatment Hierarchy & P-Score Calculations (Frequentist SUCRA)
# Output:  outputs/tables/treatment_rankings.csv
# Package: netmeta
# ==============================================================================

suppressPackageStartupMessages({
  library(netmeta)
})

cat("\n======================================================================\n")
cat(" [ANALYSIS 2/4] TREATMENT RANKING VIA P-SCORES (SUCRA ANALOGUE)\n")
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

# 3. Compute P-Scores (Frequentist Analogue to Bayesian SUCRA)
# small.values = "good" because lower Hazard Ratio represents superior overall survival
rk <- netrank(nma, small.values = "good")
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

# 4. Construct Comprehensive Treatment Ranking Table
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

# Export Ranking Table to CSV
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
output_tbl <- "outputs/tables/treatment_rankings.csv"
write.csv(df_rankings[, c("Treatment", "Rank", "Pscore_Random", "Pscore_Common", "HR_vs_Chemo_Random", "Pval_vs_Chemo")], 
          output_tbl, row.names = FALSE)

cat(sprintf(" - Successfully exported ranking table to: %s\n\n", output_tbl))
print(df_rankings[, c("Treatment", "Rank", "Pscore_Random", "HR_vs_Chemo_Random", "Pval_vs_Chemo")])
cat("\n")
