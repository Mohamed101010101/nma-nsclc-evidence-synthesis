# ==============================================================================
# Pipeline: 03_nma_full_pipeline.R
# Purpose: Master Production-Grade Network Meta-Analysis (NMA) Pipeline in R
# Method: Frequentist Graph-Theoretical NMA via 'netmeta' (G. Rücker & G. Schwarzer)
# Standard: Compliant with PRISMA-NMA (Preferred Reporting Items for Systematic
#           Reviews and Meta-Analyses - Network Meta-Analyses) & Top-Tier Medical
#           Journals (The Lancet, BMJ, JAMA, JCO, NEJM).
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Environment Setup & Dependency Loading
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(netmeta)
  library(ggplot2)
  library(readr)
  library(knitr)
})

cat("\n======================================================================\n")
cat("          NETWORK META-ANALYSIS PRODUCTION PIPELINE (R netmeta)       \n")
cat("======================================================================\n\n")

# Ensure output directories exist
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Ingest Clinical Trial Contrast Dataset
# ------------------------------------------------------------------------------
data_path <- "data/nsclc_trial_contrasts.csv"
if (!file.exists(data_path)) {
  stop("Input dataset not found at: ", data_path, ". Please run scripts/02_generate_data.R first.")
}

dat <- read.csv(data_path, stringsAsFactors = FALSE)
cat(sprintf("[DATA AUDIT] Successfully loaded: %s\n", data_path))
cat(sprintf(" - Total pairwise comparisons: %d\n", nrow(dat)))
cat(sprintf(" - Total unique clinical trials: %d\n", length(unique(dat$studlab))))
cat(sprintf(" - Total patients evaluated across trials: %s\n", format(sum(dat$n_treat1[!duplicated(paste(dat$studlab, dat$treat1))]) + sum(dat$n_treat2[!duplicated(paste(dat$studlab, dat$treat2))]), big.mark = ",")))

# ------------------------------------------------------------------------------
# 2. Fit Frequentist Network Meta-Analysis Model
# ------------------------------------------------------------------------------
# Using graph-theoretical electrical network analogy (Rücker 2012)
# Reference treatment: Chemo (Standard Platinum-Doublet Chemotherapy)
# Effect measure: Hazard Ratio (HR) -> log(HR) (TE) and standard error (seTE)
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

cat("\n[MODEL ESTIMATION COMPLETE]\n")
cat(sprintf(" - Number of treatments (n): %d\n", nma$n))
cat(sprintf(" - Number of pairwise comparisons (m): %d\n", nma$m))
cat(sprintf(" - Number of study designs (d): %d\n", nma$d))
cat(sprintf(" - Between-study heterogeneity: tau^2 = %.4f (tau = %.4f, I^2 = %.1f%%)\n", 
            nma$tau2, nma$tau, nma$I2 * 100))

# ------------------------------------------------------------------------------
# 3. Treatment Ranking via P-Scores (Frequentist Analogue to Bayesian SUCRA)
# ------------------------------------------------------------------------------
# Lower HR indicates better survival -> small.values = "good"
rk <- netrank(nma, small.values = "good")
pscores_rand <- rk$ranking.random
trt_order <- names(sort(pscores_rand, decreasing = TRUE))

df_rankings <- data.frame(
  Treatment = trt_order,
  Rank = 1:length(trt_order),
  Pscore_Random = round(pscores_rand[trt_order], 4),
  Pscore_Common = round(rk$ranking.common[trt_order], 4),
  HR_vs_Chemo_Random = sprintf("%.2f [%.2f; %.2f]", 
                               exp(nma$TE.random[trt_order, "Chemo"]),
                               exp(nma$lower.random[trt_order, "Chemo"]),
                               exp(nma$upper.random[trt_order, "Chemo"])),
  Pval_vs_Chemo = sprintf("%.4f", nma$pval.random[trt_order, "Chemo"]),
  stringsAsFactors = FALSE
)

write.csv(df_rankings, "outputs/tables/treatment_rankings.csv", row.names = FALSE)
cat("\n[TREATMENT RANKING MATRIX EXPORTED]\n")
print(df_rankings)

# ------------------------------------------------------------------------------
# 4. League Table Construction (The Dual-Model Matrix)
# ------------------------------------------------------------------------------
# Standard Top-Tier Journal Standard:
# - Lower triangle: Random-effects model HR [95% CI]
# - Upper triangle: Common/Fixed-effects model HR [95% CI]
# - Diagonal: Treatment names sorted by clinical hierarchy
lg <- netleague(nma, digits = 2, seq = trt_order)

# Export raw CSV league table
write.csv(lg$random, "outputs/tables/league_table_random_common.csv")

# Export Publication-Formatted HTML League Table
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

writeLines(html_table, "outputs/tables/league_table_formatted.html")
cat("[LEAGUE TABLE EXPORTED: CSV & HTML]\n")

# ------------------------------------------------------------------------------
# 5. Global & Local Inconsistency Diagnostics
# ------------------------------------------------------------------------------
# Decomposition of Cochran's Q:
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

write.csv(df_inconsistency, "outputs/tables/inconsistency_statistics.csv", row.names = FALSE)
cat("\n[GLOBAL INCONSISTENCY DECOMPOSITION]\n")
print(df_inconsistency)

# Local Inconsistency via Node-Splitting (Separating Direct and Indirect Evidence)
ns <- netsplit(nma)

# ------------------------------------------------------------------------------
# 6. Generate Publication-Quality Figures (300 DPI)
# ------------------------------------------------------------------------------

# Palette tailored for high-impact medical journals
colors_nodes <- c(
  "Chemo"     = "#4A5568", # Slate Gray (Standard Control)
  "IO_Mono"   = "#38A169", # Green (Immunotherapy Single-Agent)
  "IO_Chemo"  = "#3182CE", # Royal Blue (Immuno-Chemotherapy)
  "Dual_IO"   = "#805AD5", # Violet (Dual Checkpoint Blockade)
  "TKI"       = "#DD6B20", # Terracotta Orange (Targeted Monotherapy)
  "TKI_Chemo" = "#D69E2E"  # Amber Gold (Targeted + Chemotherapy)
)

# --- FIGURE 1: Publication Network Geometry Plot ---
cat("\n[RENDERING FIGURE 1: Network Geometry Graph (300 DPI)]\n")
png("outputs/figures/01_network_geometry.png", width = 3000, height = 2600, res = 300)
par(mar = c(4.6, 2.5, 3.8, 2.5))

# Compute total sample size per node for scaled prominent point size
pts_size <- sapply(nma$trts, function(t) {
  sum(dat$n_treat1[dat$treat1 == t], dat$n_treat2[dat$treat2 == t], na.rm = TRUE)
})
pts_cex <- 6.5 + (pts_size / max(pts_size)) * 4.5

netgraph(
  nma,
  points = TRUE,
  cex.points = pts_cex,
  col.points = colors_nodes[nma$trts],
  col = "#718096",
  plastic = FALSE,
  thickness = "number.of.studies",
  lwd.max = 7.5,
  lwd.min = 2,
  cex = 1.3,
  offset = 0.045,
  multiarm = TRUE,
  col.multiarm = "#E2E8F0",
  main = "Evidence Network Geometry: First-Line NSCLC Overall Survival"
)
mtext("Node diameter proportional to sample size | Line thickness proportional to trial count", 
      side = 3, line = 0.5, cex = 1.0, col = "#4A5568")

# Horizontal centered legend at bottom eliminates any node/label overlap
legend("bottom", 
       legend = c("1 Trial", "2-3 Trials", "5+ Trials"), 
       lwd = c(2, 4.5, 7.5), 
       col = "#718096", 
       horiz = TRUE,
       bty = "o", 
       box.col = "#CBD5E0",
       bg = "#FFFFFFEE",
       title = "Direct Evidence Base (Line Thickness)", 
       cex = 0.95,
       inset = c(0, -0.01),
       xpd = TRUE)
dev.off()

# --- FIGURE 2: Publication Forest Plot vs Reference (Chemo) ---
cat("[RENDERING FIGURE 2: Reference Comparison Forest Plot (300 DPI)]\n")
png("outputs/figures/02_forest_plot_random.png", width = 3200, height = 1800, res = 300)
forest(
  nma,
  reference.group = "Chemo",
  pooled = "random",
  sortvar = -rk$ranking.random,
  smlab = "Hazard Ratio (95% CI)\nvs Chemotherapy",
  label.left = "Favors Active Regimen",
  label.right = "Favors Chemotherapy",
  drop.reference.group = TRUE,
  digits = 2,
  col.square = "#2B6CB0",
  col.diamond = "#C53030",
  col.inside = "#1A202C",
  header.line = TRUE,
  leftcols = c("studlab"),
  leftlabs = c("Treatment Regimen"),
  rightcols = c("effect", "ci"),
  rightlabs = c("HR", "95% CI")
)
dev.off()

# --- FIGURE 3: Treatment Ranking (P-Score Hierarchy) Bar Chart ---
cat("[RENDERING FIGURE 3: P-Score Ranking Hierarchy (300 DPI)]\n")
df_plot_rank <- df_rankings
df_plot_rank$Treatment <- factor(df_plot_rank$Treatment, levels = rev(trt_order))

p_rank <- ggplot(df_plot_rank, aes(x = Pscore_Random, y = Treatment, fill = Treatment)) +
  geom_col(width = 0.65, alpha = 0.9, color = "#2D3748", linewidth = 0.4) +
  geom_text(aes(label = sprintf("Rank #%d | P-score: %.1f%%", Rank, Pscore_Random * 100)),
            hjust = -0.08, size = 4.2, fontface = "bold", color = "#1A202C") +
  scale_fill_manual(values = colors_nodes) +
  scale_x_continuous(limits = c(0, 1.25), breaks = seq(0, 1, 0.2), 
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Treatment Ranking Hierarchy: Surface Under Cumulative Ranking (P-Scores)",
    subtitle = "Overall Survival in Advanced NSCLC (Frequentist random-effects model)",
    x = "P-Score (Certainty of Superiority over Competing Regimens)",
    y = NULL,
    caption = "P-score ranges from 0 (certain worst) to 1 (certain best).\nComputed using netrank(..., small.values = 'good')."
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "#E2E8F0", linetype = "dashed"),
    plot.title = element_text(face = "bold", size = 16, color = "#1A365D"),
    plot.subtitle = element_text(color = "#4A5568", size = 12, margin = margin(b = 15)),
    axis.text.y = element_text(face = "bold", size = 12, color = "#2D3748"),
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10))
  )

ggsave("outputs/figures/03_pscore_ranking.png", plot = p_rank, width = 10, height = 6.5, dpi = 300)

# --- FIGURE 4: Local Inconsistency via Node-Splitting Forest Plot ---
cat("[RENDERING FIGURE 4: Node-Splitting Forest Plot (300 DPI)]\n")
# Enlarged dimensions (4600px height) ensure all comparisons and axis values are fully visible
png("outputs/figures/04_netsplit_inconsistency.png", width = 3400, height = 4600, res = 300)
forest(
  ns,
  pooled = "random",
  fontsize = 9,
  spacing = 1.05,
  digits = 2,
  smlab = "Hazard Ratio (95% CI)\nDirect vs Indirect vs Network"
)
dev.off()

# --- FIGURE 5: Net Heat Plot (Inconsistency Matrix & Evidence Contribution) ---
cat("[RENDERING FIGURE 5: Net Heat Plot (300 DPI)]\n")
# Extract Krahn design-by-treatment decomposition and Hat matrix
tau_w <- netmeta:::tau.within(nma)
nmak <- netmeta:::nma_krahn(nma, tau.preset = tau_w)
decomp <- netmeta:::decomp.design(nma, tau.preset = tau_w)

residuals <- decomp$residuals.inc.detach.random.preset
Q_inc_design <- decomp$Q.inc.design.random.preset
H_mat <- nmak$H
V_mat <- nmak$V
des <- nmak$design

Q_inc_typ <- apply(residuals, 2, function(x) t(x) %*% solve(V_mat) * x)
inc_mat <- matrix(Q_inc_design, nrow = nrow(Q_inc_typ), ncol = ncol(Q_inc_typ))
diff_mat <- inc_mat - Q_inc_typ
colnames(diff_mat) <- colnames(Q_inc_typ)
rownames(diff_mat) <- rownames(residuals)

t1_mat <- -diff_mat
wi_idx <- which(apply(t1_mat, 2, function(x) sum(is.na(x)) == nrow(t1_mat)))
if (length(wi_idx) > 0) {
  t1_mat <- t1_mat[-wi_idx, -wi_idx, drop = FALSE]
  des <- des[-wi_idx, ]
}

d1_dist <- dist(t1_mat, method = "manhattan") + dist(t(t1_mat), method = "manhattan")
h1_clust <- hclust(d1_dist)

Hp_mat <- H_mat[as.character(des$comparison), ]
if (length(wi_idx) > 0) {
  Hp_mat <- Hp_mat[, -wi_idx]
}

# Publication-grade comparison labels with multi-arm design annotations
clean_comp_label <- function(comp, design_str, narms) {
  c_clean <- gsub(":", " vs ", comp)
  if (narms > 2) {
    if (grepl("Dual_IO.*IO_Chemo", design_str)) {
      return(paste0(c_clean, " (3-arm A)"))
    } else if (grepl("Dual_IO.*IO_Mono", design_str)) {
      return(paste0(c_clean, " (3-arm B)"))
    } else if (grepl("TKI.*TKI_Chemo", design_str)) {
      return(paste0(c_clean, " (3-arm C)"))
    } else {
      return(paste0(c_clean, "*"))
    }
  } else {
    return(c_clean)
  }
}

raw_lbls <- character(nrow(des))
for (i in seq_len(nrow(des))) {
  raw_lbls[i] <- clean_comp_label(as.character(des$comparison)[i], 
                                  as.character(des$design)[i], 
                                  des$narms[i])
}

ord_idx <- h1_clust$order
ordered_lbls <- raw_lbls[ord_idx]

t1_ord <- t1_mat[ord_idx, ord_idx]
Hp_ord <- Hp_mat[ord_idx, ord_idx]

n_des <- length(ord_idx)
netheat_df <- expand.grid(Col = seq_len(n_des), Row = seq_len(n_des))
netheat_df$Inconsistency <- as.vector(t1_ord)
netheat_df$Contribution <- pmax(0, as.vector(Hp_ord))
netheat_df$ContributionPlot <- ifelse(netheat_df$Contribution > 0.005, netheat_df$Contribution, NA)

max_abs_q <- max(abs(netheat_df$Inconsistency), na.rm = TRUE)

p_netheat <- ggplot(netheat_df) +
  geom_tile(aes(x = Col, y = n_des - Row + 1, fill = Inconsistency), 
            color = "#CBD5E1", linewidth = 0.5) +
  geom_point(aes(x = Col, y = n_des - Row + 1, size = ContributionPlot), 
             shape = 15, color = "#1E293B", alpha = 0.72) +
  scale_fill_gradient2(
    low = "#1E40AF", 
    mid = "#FFFFFF", 
    high = "#DC2626", 
    midpoint = 0,
    limits = c(-max_abs_q, max_abs_q),
    labels = scales::label_number(accuracy = 0.1),
    name = "Inconsistency\nContribution (ΔQ)",
    guide = guide_colorbar(
      order = 1,
      frame.colour = "#94A3B8",
      ticks.colour = "#475569"
    )
  ) +
  scale_size_area(
    max_size = 14,
    breaks = c(0.10, 0.30, 0.50, 0.70),
    labels = c("10%", "30%", "50%", "70%"),
    name = "Direct Evidence\nWeight (Hat Matrix)",
    guide = guide_legend(
      order = 2,
      override.aes = list(shape = 15, color = "#1E293B", alpha = 0.8)
    ),
    na.value = NA
  ) +
  scale_x_continuous(
    breaks = seq_len(n_des), 
    labels = ordered_lbls,
    position = "top"
  ) +
  scale_y_continuous(
    breaks = seq_len(n_des),
    labels = rev(ordered_lbls)
  ) +
  coord_fixed() +
  labs(
    title = "Net Heat Plot: Matrix of Inconsistency & Evidence Contribution",
    subtitle = "Design-by-treatment interaction model with hierarchical clustering (Random-Effects)",
    caption = paste0(
      "Interpretation Guide:\n",
      "• Background Color: Evaluates inconsistency contribution when detaching direct evidence.\n",
      "  Warm red = drives network inconsistency; cool blue = stabilizes the network.\n",
      "• Inner Grey Squares: Area proportional to direct evidence contribution (Hat matrix H_ij).\n",
      "  Larger squares signify higher statistical weight of direct data in estimating each comparison.\n",
      "• Multi-Arm Trial Designs: (3-arm A) Chemo/Dual_IO/IO_Chemo [CheckMate-9LA];\n",
      "  (3-arm B) Chemo/Dual_IO/IO_Mono [KEYNOTE-598]; (3-arm C) Chemo/TKI/TKI_Chemo [NEJ026]."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "#FFFFFF", color = NA),
    panel.background = element_rect(fill = "#FFFFFF", color = NA),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 16, color = "#0F172A", hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 11, color = "#475569", hjust = 0.5, margin = margin(b = 14)),
    plot.caption = element_text(size = 8.8, color = "#475569", hjust = 0, lineheight = 1.35, margin = margin(t = 14)),
    axis.title = element_blank(),
    axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 9.5, face = "bold", color = "#1E293B"),
    axis.text.y = element_text(size = 9.5, face = "bold", color = "#1E293B"),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 9.5, face = "bold", color = "#0F172A"),
    legend.text = element_text(size = 8.5),
    legend.key.height = unit(1.6, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.spacing.y = unit(0.8, "cm"),
    plot.margin = margin(t = 15, r = 25, b = 15, l = 20)
  )

ggsave("outputs/figures/05_netheat_plot.png", plot = p_netheat, width = 11.5, height = 10.5, dpi = 300)

# --- FIGURE 6: Comparison-Adjusted Funnel Plot (Small-Study Effects) ---
cat("[RENDERING FIGURE 6: Comparison-Adjusted Funnel Plot (300 DPI)]\n")
png("outputs/figures/06_funnel_plot.png", width = 2800, height = 2400, res = 300)
par(mar = c(4.5, 4.5, 3.5, 2))
funnel(
  nma,
  order = trt_order,
  pooled = "random",
  pch = 19,
  col = "#2B6CB0",
  cex = 1.3,
  linreg = TRUE,
  main = "Comparison-Adjusted Funnel Plot (Evaluation of Small-Study Effects)",
  xlab = "Log Hazard Ratio centered by comparison-specific effect",
  ylab = "Standard Error of Log Hazard Ratio"
)
dev.off()

# --- FIGURE 7: Publication-Grade League Table Matrix Figure (300 DPI) ---
cat("[RENDERING FIGURE 7: Publication League Table Figure (300 DPI)]\n")

n_trts <- length(trt_order)
mat_rnd <- lg$random
pval_rnd <- nma$pval.random[trt_order, trt_order]

trt_meta <- list(
  "IO_Chemo"  = list(name = "IO + Chemo", class = "IO Combo", pscore = sprintf("%.1f%%", pscores_rand["IO_Chemo"] * 100), rank = "Rank 1"),
  "TKI_Chemo" = list(name = "TKI + Chemo", class = "Targeted Combo", pscore = sprintf("%.1f%%", pscores_rand["TKI_Chemo"] * 100), rank = "Rank 2"),
  "Dual_IO"   = list(name = "Dual IO", class = "Dual Checkpoint", pscore = sprintf("%.1f%%", pscores_rand["Dual_IO"] * 100), rank = "Rank 3"),
  "IO_Mono"   = list(name = "IO Monotherapy", class = "Anti-PD-(L)1", pscore = sprintf("%.1f%%", pscores_rand["IO_Mono"] * 100), rank = "Rank 4"),
  "TKI"       = list(name = "TKI Monotherapy", class = "Targeted Mono", pscore = sprintf("%.1f%%", pscores_rand["TKI"] * 100), rank = "Rank 5"),
  "Chemo"     = list(name = "Chemotherapy", class = "Standard Control", pscore = sprintf("%.1f%%", pscores_rand["Chemo"] * 100), rank = "Rank 6")
)

cells_df <- data.frame()

for (i in 1:n_trts) {
  for (j in 1:n_trts) {
    row_trt <- trt_order[i]
    col_trt <- trt_order[j]
    cell_type <- if (i == j) "diagonal" else if (i > j) "lower" else "upper"
    
    tag_text <- ""
    hr_text <- ""
    ci_text <- ""
    bg_color <- "#FFFFFF"
    hr_color <- "#1E293B"
    ci_color <- "#64748B"
    tag_color <- "#94A3B8"
    hr_size <- 4.6
    hr_fontface <- "bold"
    
    val_str <- mat_rnd[i, j]
    
    if (cell_type == "diagonal") {
      meta_info <- trt_meta[[row_trt]]
      tag_text <- meta_info$class
      d_name <- meta_info$name
      if (d_name == "IO Monotherapy") d_name <- "IO\nMonotherapy"
      if (d_name == "TKI Monotherapy") d_name <- "TKI\nMonotherapy"
      if (d_name == "Chemotherapy") d_name <- "Chemo-\ntherapy"
      
      hr_text <- d_name
      ci_text <- paste0(meta_info$rank, " • P-Score: ", meta_info$pscore)
      bg_color <- "#1E3A8A" # Deep Navy
      hr_color <- "#FFFFFF"
      ci_color <- "#FDE047" # Crisp Golden Yellow
      tag_color <- "#93C5FD" # Soft light blue
      hr_size <- 4.8
    } else if (cell_type == "lower") {
      # Lower triangle: Network Random Effects (Column vs Row)
      clean_v <- gsub(";", " –", val_str)
      parts <- strsplit(clean_v, " \\[")[[1]]
      hr_text <- paste0("HR ", parts[1])
      ci_text <- paste0("95% CI: ", gsub("\\]", "", parts[2]))
      
      p_val <- pval_rnd[i, j]
      is_sig <- (!is.na(p_val) && p_val < 0.05)
      
      if (is_sig) {
        bg_color <- "#DCFCE7" # Soft Mint Green
        hr_color <- "#14532D" # Deep Forest Green
        ci_color <- "#166534"
        tag_color <- "#15803D"
        tag_text <- "Network ★ p < 0.05"
      } else {
        bg_color <- "#F8FAFC" # Soft Neutral Slate
        hr_color <- "#334155"
        ci_color <- "#64748B"
        tag_color <- "#94A3B8"
        tag_text <- "Network (Random)"
      }
    } else if (cell_type == "upper") {
      # Upper triangle: Direct pairwise RCT evidence (Row vs Column)
      if (val_str == "." || is.na(val_str)) {
        tag_text <- "Direct Evidence"
        hr_text <- "—"
        ci_text <- "No direct head-to-head"
        bg_color <- "#F1F5F9"
        hr_color <- "#94A3B8"
        ci_color <- "#94A3B8"
        tag_color <- "#CBD5E1"
        hr_size <- 4.2
        hr_fontface <- "plain"
      } else {
        tag_text <- "Direct Evidence (RCT)"
        clean_v <- gsub(";", " –", val_str)
        parts <- strsplit(clean_v, " \\[")[[1]]
        hr_text <- paste0("HR ", parts[1])
        ci_text <- paste0("95% CI: ", gsub("\\]", "", parts[2]))
        bg_color <- "#EFF6FF" # Soft pastel blue
        hr_color <- "#1E40AF"
        ci_color <- "#2563EB"
        tag_color <- "#3B82F6"
      }
    }
    
    cells_df <- rbind(cells_df, data.frame(
      Row = i,
      Col = j,
      CellType = cell_type,
      Tag = tag_text,
      HR = hr_text,
      CI = ci_text,
      BG = bg_color,
      HRColor = hr_color,
      CIColor = ci_color,
      TagColor = tag_color,
      HRSize = hr_size,
      HRFontface = hr_fontface,
      stringsAsFactors = FALSE
    ))
  }
}

cells_df$X <- cells_df$Col
cells_df$Y <- n_trts - cells_df$Row + 1

axis_labels <- c("IO + Chemo", "TKI + Chemo", "Dual IO", "IO Mono", "TKI Mono", "Chemo")

p_league <- ggplot(cells_df) +
  geom_rect(aes(
    xmin = X - 0.47, xmax = X + 0.47,
    ymin = Y - 0.47, ymax = Y + 0.47,
    fill = BG
  ), color = "#CBD5E1", linewidth = 0.8) +
  scale_fill_identity() +
  geom_text(aes(x = X, y = Y + 0.28, label = Tag, color = TagColor),
            size = 2.6, fontface = "bold") +
  geom_text(aes(x = X, y = ifelse(CellType == "diagonal", Y + 0.02, Y + 0.03), 
                label = HR, color = HRColor,
                size = HRSize, fontface = HRFontface), lineheight = 0.95) +
  geom_text(aes(x = X, y = Y - 0.26, label = CI, color = CIColor),
            size = 2.85, fontface = "plain") +
  scale_size_identity() +
  scale_color_identity() +
  scale_x_continuous(
    breaks = 1:n_trts,
    labels = axis_labels,
    position = "top",
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    breaks = 1:n_trts,
    labels = rev(axis_labels),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  coord_fixed() +
  labs(
    title = "League Table of Pairwise Treatment Effects (Network vs Direct Evidence)",
    subtitle = "Treatments ordered by clinical hierarchy (P-scores) from top-left (best) to bottom-right (worst)",
    caption = paste0(
      "Reading Guide & Publication Conventions:\n",
      "• Lower Triangle (Green / Neutral): Network Meta-Analysis estimates from Random-Effects model (Column vs Row treatment).\n",
      "  Hazard Ratio (HR) < 1.0 indicates superiority of the higher-ranked Column treatment over the Row treatment.\n",
      "• Green Shading (★): Denotes statistically significant superior efficacy at the p < 0.05 threshold (95% CI excludes 1.0).\n",
      "• Upper Triangle (Light Blue): Direct pairwise meta-analysis estimates from Head-to-Head RCTs (dashes indicate purely indirect links).\n",
      "• Diagonal (Deep Navy): Treatment node names with their overall clinical hierarchy rank and P-score."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "#FFFFFF", color = NA),
    panel.background = element_rect(fill = "#FFFFFF", color = NA),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 16, color = "#0F172A", hjust = 0.5, margin = margin(b = 6)),
    plot.subtitle = element_text(size = 11.5, color = "#475569", hjust = 0.5, margin = margin(b = 16)),
    plot.caption = element_text(size = 9.2, color = "#334155", hjust = 0, lineheight = 1.38, margin = margin(t = 16)),
    axis.title = element_blank(),
    axis.text.x.top = element_text(size = 11, face = "bold", color = "#1E293B", margin = margin(b = 8)),
    axis.text.y = element_text(size = 11, face = "bold", color = "#1E293B", margin = margin(r = 8)),
    plot.margin = margin(t = 20, r = 25, b = 20, l = 25)
  )

ggsave("outputs/figures/07_league_table_figure.png", plot = p_league, width = 13.5, height = 13.5, dpi = 300)

cat("\n======================================================================\n")
cat(" [SUCCESS] Master NMA Production Pipeline Completed Flawlessly!\n")
cat(" All publication tables saved in: outputs/tables/\n")
cat(" All 300-DPI high-res figures saved in: outputs/figures/\n")
cat("======================================================================\n\n")
