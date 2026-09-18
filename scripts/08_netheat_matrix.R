# ==============================================================================
# Script: 08_netheat_matrix.R
# Purpose: Net Heat Plot (Matrix of Inconsistency & Direct Evidence Contribution)
# Output:  outputs/figures/05_netheat_plot.png (300 DPI Publication Figure)
# Package: netmeta & ggplot2 (Krahn Design-by-Treatment Interaction Model)
# ==============================================================================

suppressPackageStartupMessages({
  library(netmeta)
  library(ggplot2)
  library(scales)
})

cat("\n======================================================================\n")
cat(" [ANALYSIS 6/7] NET HEAT MATRIX: INCONSISTENCY & HAT MATRIX WEIGHTS\n")
cat("======================================================================\n")

# 1. Load Clinical Trial Contrast Data
data_path <- "data/nsclc_trial_contrasts.csv"
if (!file.exists(data_path)) {
  stop(sprintf("Data file not found at: %s. Please run scripts/02_generate_data.R first.", data_path))
}
dat <- read.csv(data_path, stringsAsFactors = FALSE)
cat(sprintf(" - Loaded contrast dataset: %d comparisons across %d trials\n",
            nrow(dat), length(unique(dat$studlab))))

# 2. Fit Frequentist Graph-Theoretical Model
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

# 3. Extract Krahn Design-by-Treatment Decomposition and Hat Matrix
cat(" - Computing Krahn decomposition and Hat matrix weights ...\n")
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

# 4. Format Publication-Grade Comparison Labels with Multi-Arm Annotations
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

# 5. Render Publication Net Heat Plot (300 DPI)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
output_fig <- "outputs/figures/05_netheat_plot.png"
cat(sprintf(" - Rendering Figure 5 to: %s ...\n", output_fig))

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

ggsave(output_fig, plot = p_netheat, width = 11.5, height = 10.5, dpi = 300)

cat(sprintf(" [SUCCESS] Net Heat plot saved cleanly: %s\n\n", output_fig))
