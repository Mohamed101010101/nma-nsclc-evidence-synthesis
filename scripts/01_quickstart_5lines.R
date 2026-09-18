# ==============================================================================
# Script: 01_quickstart_5lines.R
# Purpose: The 5-Line Pedagogical Quickstart to Network Meta-Analysis in R
# "Why start here? Because with just 5 lines of R, you get a full network and
#  league table, intuitively grasping indirect comparisons before Bayesian math."
# ==============================================================================

# Line 1: Load the frequentist network meta-analysis engine
library(netmeta)

# Line 2: Read the contrast-level clinical trial dataset (ln(HR) and seTE)
dat <- read.csv("data/nsclc_trial_contrasts.csv")

# Line 3: Fit the full Network Meta-Analysis (Random & Common Effects)
nma <- netmeta(TE, seTE, treat1, treat2, studlab, data = dat, sm = "HR", ref = "Chemo")

# Line 4: Visualize the complete Network Geometry
netgraph(nma, points = TRUE, cex = 1.5, col = "#1f77b4", plastic = FALSE, thickness = "number.of.studies")

# Line 5: Print the League Table of all pairwise comparisons (Random-Effects)
print(netleague(nma, digits = 2)$random)
