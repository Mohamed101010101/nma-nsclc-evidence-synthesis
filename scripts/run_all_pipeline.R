# ==============================================================================
# Script: run_all_pipeline.R
# Purpose: Master Orchestrator for All NMA Analyses & Publication Figure Designs
# Architecture:
#   - Phase 1: Statistical Analyses & Table Generation (scripts/analyses/)
#   - Phase 2: Publication Figure Designs (scripts/designs/)
# Execution: Rscript scripts/run_all_pipeline.R
# ==============================================================================

cat("\n==============================================================================\n")
cat("          MASTER NETWORK META-ANALYSIS PRODUCTION PIPELINE                    \n")
cat("          Advanced Non-Small Cell Lung Cancer (NSCLC) Evidence Synthesis      \n")
cat("==============================================================================\n")

start_time <- Sys.time()

# 0. Ensure contrast dataset exists
if (!file.exists("data/nsclc_trial_contrasts.csv")) {
  stop("Error: Contrast dataset 'data/nsclc_trial_contrasts.csv' not found. Please ensure the data file exists before running the pipeline.")
}

# ------------------------------------------------------------------------------
# PHASE 1: STATISTICAL ANALYSES & TABLE GENERATION
# ------------------------------------------------------------------------------
cat("\n==============================================================================\n")
cat(" >>> PHASE 1: STATISTICAL ANALYSES & SUMMARY TABLES (scripts/analyses/)        \n")
cat("==============================================================================\n")

analysis_modules <- list(
  list(file = "scripts/analyses/01_fit_nma_model.R",      name = "Model Estimation & Heterogeneity Assessment"),
  list(file = "scripts/analyses/02_treatment_rankings.R", name = "Treatment Hierarchy & P-Scores (outputs/tables/treatment_rankings.csv)"),
  list(file = "scripts/analyses/03_league_table.R",       name = "Dual-Model League Table (CSV Matrix)"),
  list(file = "scripts/analyses/04_inconsistency_tests.R",name = "Global Q Decomposition & Inconsistency Table"),
  list(file = "scripts/analyses/05_league_table_html.R",  name = "Formatted Interactive HTML League Table")
)

for (idx in seq_along(analysis_modules)) {
  mod <- analysis_modules[[idx]]
  cat(sprintf("\n[ANALYSIS %d/%d]: %s (%s) ...\n", idx, length(analysis_modules), mod$name, mod$file))
  t0 <- Sys.time()
  
  env <- new.env(parent = globalenv())
  source(mod$file, local = env)
  
  t1 <- Sys.time()
  cat(sprintf(">>> Analysis %d finished in %.2f seconds.\n", idx, as.numeric(difftime(t1, t0, units = "secs"))))
}

# ------------------------------------------------------------------------------
# PHASE 2: PUBLICATION GRAPHIC DESIGNS (FIGURES 1 TO 7)
# ------------------------------------------------------------------------------
cat("\n==============================================================================\n")
cat(" >>> PHASE 2: PUBLICATION GRAPHIC DESIGNS (300 DPI) (scripts/designs/)        \n")
cat("==============================================================================\n")

design_modules <- list(
  list(file = "scripts/designs/fig01_network_geometry.R",     name = "Figure 1: Evidence Network Geometry (Topology)"),
  list(file = "scripts/designs/fig02_forest_plot.R",         name = "Figure 2: Reference Forest Plot vs Chemotherapy"),
  list(file = "scripts/designs/fig03_pscore_ranking.R",      name = "Figure 3: P-Score Treatment Ranking Hierarchy"),
  list(file = "scripts/designs/fig04_netsplit_inconsistency.R", name = "Figure 4: Node-Splitting Local Inconsistency Forest Plot"),
  list(file = "scripts/designs/fig05_netheat_plot.R",        name = "Figure 5: Net Heat Inconsistency Matrix Plot"),
  list(file = "scripts/designs/fig06_funnel_plot.R",         name = "Figure 6: Comparison-Adjusted Funnel Plot"),
  list(file = "scripts/designs/fig07_league_table_matrix.R", name = "Figure 7: Publication League Table Graphic Matrix")
)

for (idx in seq_along(design_modules)) {
  mod <- design_modules[[idx]]
  cat(sprintf("\n[DESIGN %d/%d]: %s (%s) ...\n", idx, length(design_modules), mod$name, mod$file))
  t0 <- Sys.time()
  
  env <- new.env(parent = globalenv())
  source(mod$file, local = env)
  
  t1 <- Sys.time()
  cat(sprintf(">>> Design %d finished in %.2f seconds.\n", idx, as.numeric(difftime(t1, t0, units = "secs"))))
}

total_elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat("\n==============================================================================\n")
cat(sprintf(" [PIPELINE SUCCESS] All %d analyses and %d figure designs generated in %.2f seconds!\n", 
            length(analysis_modules), length(design_modules), total_elapsed))
cat(" Outputs Directory Audit:\n")
cat(sprintf("  - Figures: outputs/figures/ (%d publication-grade 300 DPI figures, Fig 01 to %02d)\n",
            length(design_modules), length(design_modules)))
cat("  - Tables:  outputs/tables/  (3 CSVs + 1 Formatted Interactive HTML Table)\n")
cat("==============================================================================\n\n")
