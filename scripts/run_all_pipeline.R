# ==============================================================================
# Script: run_all_pipeline.R
# Purpose: Master Orchestrator for All Modular Network Meta-Analysis Scripts
# Execution: Rscript scripts/run_all_pipeline.R
# ==============================================================================

cat("\n==============================================================================\n")
cat("          STARTING COMPREHENSIVE NETWORK META-ANALYSIS PIPELINE               \n")
cat("          Advanced Non-Small Cell Lung Cancer (NSCLC) Evidence Synthesis      \n")
cat("==============================================================================\n")

start_time <- Sys.time()

# Verify that dataset exists, or generate it automatically
if (!file.exists("data/nsclc_trial_contrasts.csv")) {
  cat("\n[STEP 0] Dataset not detected. Generating synthetic trial contrast data ...\n")
  source("scripts/02_generate_data.R", local = new.env())
}

# Define modular analysis modules
modules <- list(
  list(file = "scripts/03_network_geometry.R",     name = "Network Geometry Graph (Figure 1)"),
  list(file = "scripts/04_forest_plot.R",         name = "Reference Forest Plot vs Chemo (Figure 2)"),
  list(file = "scripts/05_treatment_rankings.R",  name = "P-Score Rankings & Hierarchy (Figure 3 + CSV)"),
  list(file = "scripts/06_league_table.R",        name = "Dual-Model League Table (Figure 7 + CSV + HTML)"),
  list(file = "scripts/07_inconsistency_netsplit.R", name = "Global Q Decomposition & Node-Splitting (Figure 4 + CSV)"),
  list(file = "scripts/08_netheat_matrix.R",      name = "Net Heat Inconsistency Matrix (Figure 5)"),
  list(file = "scripts/09_funnel_plot.R",         name = "Comparison-Adjusted Funnel Plot (Figure 6)")
)

total_steps <- length(modules)

for (idx in seq_along(modules)) {
  mod <- modules[[idx]]
  cat(sprintf("\n>>> Executing [%d/%d]: %s (%s) ...\n", idx, total_steps, mod$name, mod$file))
  t0 <- Sys.time()
  
  # Run module in isolated environment to avoid cross-script contamination
  env <- new.env(parent = globalenv())
  source(mod$file, local = env)
  
  t1 <- Sys.time()
  cat(sprintf(">>> Completed [%d/%d] in %.2f seconds.\n", idx, total_steps, as.numeric(difftime(t1, t0, units = "secs"))))
}

total_elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat("\n==============================================================================\n")
cat(sprintf(" [PIPELINE SUCCESS] All %d analyses and publication figures generated in %.2f seconds!\n", 
            total_steps, total_elapsed))
cat(" Outputs Directory Audit:\n")
cat("  - Figures: outputs/figures/ (7 publication-grade 300 DPI figures)\n")
cat("  - Tables:  outputs/tables/  (3 CSVs + 1 Formatted Interactive HTML)\n")
cat("==============================================================================\n\n")
