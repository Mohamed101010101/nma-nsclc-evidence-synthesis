# ==============================================================================
# Script: scripts/analyses/03_league_table.R
# Purpose: Pairwise Dual-Model League Table Generation (CSV & Publication HTML)
# Outputs: outputs/tables/league_table_random_common.csv
#          outputs/tables/league_table_formatted.html
# Package: netmeta
# ==============================================================================

suppressPackageStartupMessages({
  library(netmeta)
})

cat("\n======================================================================\n")
cat(" [ANALYSIS 3/4] LEAGUE TABLE GENERATION (CSV & FORMATTED HTML)\n")
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
cat(sprintf(" - Exported raw league matrix CSV to: %s\n", output_csv))

# 5. Export Publication-Formatted HTML League Table
output_html <- "outputs/tables/league_table_formatted.html"
html_table <- paste0(
  "<div style='font-family: Arial, sans-serif; margin: 20px 0;'>\n",
  "<h3 style='color: #1a365d; text-align: center;'>Table: League Table of Pairwise Treatment Comparisons (Hazard Ratios [95% CI])</h3>\n",
  "<p style='text-align: center; color: #4a5568; font-size: 0.9em;'>Treatments ordered by hierarchy (P-scores) from top-left (best) to bottom-right (worst).<br>",
  "<b>Lower Triangle:</b> Random-Effects Model | <b>Upper Triangle:</b> Common-Effects Model | <b>Bold:</b> Significant (p < 0.05)</p>\n",
  "<table style='border-collapse: collapse; margin: 0 auto; width: 95%; box-shadow: 0 4px 6px rgba(0,0,0,0.1); border-radius: 8px; overflow: hidden;'>\n",
  "  <thead>\n    <tr style='background-color: #2b6cb0; color: white; text-align: center; font-weight: bold;'>\n",
  "      <th style='padding: 12px; border: 1px solid #cbd5e0;'>Treatment</th>\n",
  paste0(sprintf("      <th style='padding: 12px; border: 1px solid #cbd5e0;'>%s</th>\n", trt_order), collapse = ""),
  "    </tr>\n  </thead>\n  <tbody>\n"
)

mat <- lg$random
for (i in 1:nrow(mat)) {
  row_html <- sprintf("    <tr style='background-color: %s; text-align: center;'>\n", ifelse(i %% 2 == 0, "#f7fafc", "#ffffff"))
  row_html <- paste0(row_html, sprintf("      <td style='padding: 10px; font-weight: bold; background-color: #edf2f7; border: 1px solid #cbd5e0;'>%s</td>\n", rownames(mat)[i]))
  for (j in 1:ncol(mat)) {
    val <- mat[i, j]
    is_diag <- (i == j)
    cell_style <- if (is_diag) {
      "padding: 10px; font-weight: bold; background-color: #bee3f8; color: #2b6cb0; border: 1px solid #cbd5e0;"
    } else {
      "padding: 10px; border: 1px solid #cbd5e0; font-size: 0.95em;"
    }
    row_html <- paste0(row_html, sprintf("      <td style='%s'>%s</td>\n", cell_style, val))
  }
  row_html <- paste0(row_html, "    </tr>\n")
  html_table <- paste0(html_table, row_html)
}

html_table <- paste0(
  html_table,
  "  </tbody>\n</table>\n",
  "<p style='font-size: 0.85em; color: #718096; text-align: center; margin-top: 8px;'>",
  "HR < 1 favors column-defining treatment in lower triangle, and row-defining treatment in upper triangle.",
  "</p>\n</div>"
)

writeLines(html_table, output_html)
cat(sprintf(" - Exported formatted HTML league table to: %s\n\n", output_html))
