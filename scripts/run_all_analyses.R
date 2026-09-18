# ==============================================================================
# Script: run_all_analyses.R
# Purpose: Dedicated Runner for All 5 NMA Statistical Analyses & Tables
# Execution: Rscript scripts/run_all_analyses.R
# ==============================================================================

cat("\n==============================================================================\n")
cat("          RUNNING ALL STATISTICAL ANALYSES (scripts/analyses/)                \n")
cat("          Advanced Non-Small Cell Lung Cancer (NSCLC) Evidence Synthesis      \n")
cat("==============================================================================\n")

start_time <- Sys.time()

# 0. Ensure contrast dataset exists
if (!file.exists("data/nsclc_trial_contrasts.csv")) {
  cat("\n[STEP 0] Generating synthetic trial contrast data ...\n")
  source("scripts/02_generate_data.R", local = new.env())
}

analysis_modules <- list(
  list(file = "scripts/analyses/01_fit_nma_model.R",       name = "1. Model Estimation & Heterogeneity Assessment"),
  list(file = "scripts/analyses/02_treatment_rankings.R",  name = "2. Treatment Hierarchy & P-Scores"),
  list(file = "scripts/analyses/03_league_table.R",        name = "3. Dual-Model League Table (CSV Matrix)"),
  list(file = "scripts/analyses/04_inconsistency_tests.R", name = "4. Global Q Decomposition & Inconsistency Tests"),
  list(file = "scripts/analyses/05_league_table_html.R",   name = "5. Formatted Interactive HTML League Table")
)

for (idx in seq_along(analysis_modules)) {
  mod <- analysis_modules[[idx]]
  cat(sprintf("\n------------------------------------------------------------------------------\n"))
  cat(sprintf(" >>> EXECUTING: %s (%s) ...\n", mod$name, mod$file))
  cat(sprintf("------------------------------------------------------------------------------\n"))
  t0 <- Sys.time()
  
  env <- new.env(parent = globalenv())
  source(mod$file, local = env)
  
  t1 <- Sys.time()
  cat(sprintf(">>> Analysis %d finished in %.2f seconds.\n", idx, as.numeric(difftime(t1, t0, units = "secs"))))
}

total_duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat("\n==============================================================================\n")
cat(sprintf(" [SUCCESS] All 5 statistical analyses completed in %.2f seconds!\n", total_duration))
cat(" Generated Tables & Models:\n")
cat("  - Models: outputs/models/nma_model.rds, nma_rankings.rds\n")
cat("  - Table:  outputs/tables/treatment_rankings.csv\n")
cat("  - Table:  outputs/tables/league_table_random_common.csv\n")
cat("  - Table:  outputs/tables/inconsistency_statistics.csv\n")
cat("  - HTML:   outputs/tables/league_table_formatted.html\n")
cat("==============================================================================\n\n")
