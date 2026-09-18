# ==============================================================================
# Script: scripts/analyses/06_rob2_analysis.R
# Purpose: Cochrane Risk of Bias 2 (RoB 2) Statistical Synthesis & Summary Metrics
# Standard: Cochrane Handbook for Systematic Reviews of Interventions (RoB 2 Tool)
# Inputs:  data/nsclc_rob2_assessments.csv
# Outputs: outputs/models/rob2_data.rds
#          outputs/tables/rob2_study_assessments.csv
#          outputs/tables/rob2_domain_summary.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
})

cat("\n======================================================================\n")
cat(" [ANALYSIS 6/7] COCHRANE RISK OF BIAS 2 (RoB 2) STATISTICAL SYNTHESIS\n")
cat("======================================================================\n")

# 1. Load RoB 2 Assessment Data
data_path <- "data/nsclc_rob2_assessments.csv"
if (!file.exists(data_path)) {
  stop(sprintf("RoB 2 assessment dataset not found at: %s", data_path))
}

rob_df <- read.csv(data_path, stringsAsFactors = FALSE)
cat(sprintf(" - Loaded RoB 2 assessments: %d clinical trials\n", nrow(rob_df)))

domains <- c("D1", "D2", "D3", "D4", "D5", "Overall")
domain_labels <- c(
  "D1"      = "D1: Randomisation process",
  "D2"      = "D2: Deviations from intended interventions",
  "D3"      = "D3: Missing outcome data",
  "D4"      = "D4: Measurement of the outcome",
  "D5"      = "D5: Selection of the reported result",
  "Overall" = "Overall Risk of Bias"
)

levels_rob <- c("Low risk", "Some concerns", "High risk")

# 2. Compute Domain Summary Frequencies & Percentages
summary_list <- list()

for (d in domains) {
  counts <- table(factor(rob_df[[d]], levels = levels_rob))
  pcts   <- round(prop.table(counts) * 100, 1)
  
  summary_list[[d]] <- data.frame(
    Domain_Code = d,
    Domain_Name = domain_labels[d],
    Low_Risk_N  = as.integer(counts["Low risk"]),
    Low_Risk_Pct = as.numeric(pcts["Low risk"]),
    Some_Concerns_N = as.integer(counts["Some concerns"]),
    Some_Concerns_Pct = as.numeric(pcts["Some concerns"]),
    High_Risk_N = as.integer(counts["High risk"]),
    High_Risk_Pct = as.numeric(pcts["High risk"]),
    Total_Studies = nrow(rob_df),
    stringsAsFactors = FALSE
  )
}

df_summary <- do.call(rbind, summary_list)
rownames(df_summary) <- NULL

# 3. Export Summary Tables & Cache RDS Object
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/models", recursive = TRUE, showWarnings = FALSE)

study_csv   <- "outputs/tables/rob2_study_assessments.csv"
summary_csv <- "outputs/tables/rob2_domain_summary.csv"
rds_file    <- "outputs/models/rob2_data.rds"

write.csv(rob_df, study_csv, row.names = FALSE)
write.csv(df_summary, summary_csv, row.names = FALSE)

rob_bundle <- list(
  raw_data       = rob_df,
  domain_summary = df_summary,
  domain_labels  = domain_labels,
  levels         = levels_rob,
  n_studies      = nrow(rob_df)
)

saveRDS(rob_bundle, rds_file)

cat(sprintf(" - Saved study-level RoB 2 table: %s\n", study_csv))
cat(sprintf(" - Saved domain percentage summary: %s\n", summary_csv))
cat(sprintf(" - Serialized RoB 2 data bundle: %s\n", rds_file))

cat("\n [RoB 2 DOMAIN DISTRIBUTION SUMMARY]\n")
for (i in 1:nrow(df_summary)) {
  cat(sprintf(" - %-42s : Low %4.1f%% | Some Concerns %4.1f%% | High %4.1f%%\n",
              df_summary$Domain_Name[i],
              df_summary$Low_Risk_Pct[i],
              df_summary$Some_Concerns_Pct[i],
              df_summary$High_Risk_Pct[i]))
}

cat("\n [SUCCESS] Cochrane RoB 2 statistical synthesis completed successfully.\n\n")
