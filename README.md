# 🩺 Clinical Evidence Synthesis & Network Meta-Analysis Engine in R (`netmeta`)
### *Frequentist Graph-Theoretical Synthesis, Transitivity Diagnostics & Top-Tier Publication Pipeline*

[![R Version](https://img.shields.io/badge/R-v4.6.1-276DC3.svg?logo=R&logoColor=white)](https://www.r-project.org/)
[![Package: netmeta](https://img.shields.io/badge/netmeta-v3.6--1-blue.svg)](https://cran.r-project.org/package=netmeta)
[![Methodology](https://img.shields.io/badge/Methodology-10%20Advanced%20Engines-darkgreen.svg)](#-core-biostatistical--methodological-competencies)
[![Guideline: PRISMA-NMA](https://img.shields.io/badge/PRISMA--NMA-100%25%20Compliant-success.svg)](http://www.prisma-statement.org/Extensions/NetworkMetaAnalysis)
[![Standard](https://img.shields.io/badge/Journal%20Standard-Lancet%20%7C%20NEJM%20%7C%20BMJ%20%7C%20JAMA-purple.svg)](#-publication-gallery-300-dpi-visual-exhibits)
[![Evidence Base](https://img.shields.io/badge/Evidence%20Base-24%20RCTs%20%7C%2015%2C753%20Pts-informational.svg)](#-clinical-research-scenario-advanced-nsclc)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## 📑 Table of Contents
- [Executive Overview & Pedagogical Rationale](#-executive-overview--pedagogical-rationale)
- [Core Biostatistical & Methodological Competencies](#-core-biostatistical--methodological-competencies)
- [Clinical Research Scenario (Advanced NSCLC)](#-clinical-research-scenario-advanced-nsclc)
- [Mathematical & Biostatistical Framework](#-mathematical--biostatistical-framework)
  - [1. Survival Contrast Representation ($\ln(\text{HR})$, $\text{SE}$)](#1-survival-contrast-representation-lnhr-textse)
  - [2. Bucher's Principle of Indirect Comparison](#2-buchers-principle-of-indirect-comparison)
  - [3. Electrical Network Analogy & Laplacian Matrices (Rücker 2012)](#3-electrical-network-analogy--laplacian-matrices-rücker-2012)
  - [4. Frequentist P-Scores & 10,000 Monte Carlo SUCRA](#4-frequentist-p-scores--10000-monte-carlo-sucra)
  - [5. Component Network Meta-Analysis (CNMA Synergy Testing)](#5-component-network-meta-analysis-cnma-synergy-testing)
  - [6. Bi-dimensional Benefit-Risk Trade-Off Synthesis](#6-bi-dimensional-benefit-risk-trade-off-synthesis)
- [Publication Gallery (300 DPI Visual Exhibits: 12 Figures)](#-publication-gallery-300-dpi-visual-exhibits)
  - [Figure 01: Evidence Network Geometry](#figure-01-evidence-network-geometry)
  - [Figure 02: Reference Forest Plot vs Chemotherapy](#figure-02-reference-forest-plot-vs-chemotherapy)
  - [Figure 03: P-Score Treatment Ranking Hierarchy](#figure-03-p-score-treatment-ranking-hierarchy)
  - [Figure 04: Node-Splitting Local Inconsistency (`netsplit`)](#figure-04-node-splitting-local-inconsistency-netsplit)
  - [Figure 05: Net Heat Inconsistency Matrix Plot](#figure-05-net-heat-inconsistency-matrix-plot)
  - [Figure 06: Comparison-Adjusted Funnel Plot](#figure-06-comparison-adjusted-funnel-plot)
  - [Figure 07: Dual-Model League Table Matrix](#figure-07-dual-model-league-table-matrix)
  - [Figure 08: Leave-One-Out (LOO) Influence Forest Plot](#figure-08-leave-one-out-loo-influence-forest-plot)
  - [Figure 09: Component NMA (CNMA) Incremental Effects Forest](#figure-09-component-nma-cnma-incremental-effects-forest)
  - [Figure 10: Probabilistic Hierarchy & Cumulative Rankograms](#figure-10-probabilistic-hierarchy--cumulative-rankograms)
  - [Figure 11: Bi-dimensional Benefit-Risk Trade-Off Matrix](#figure-11-bi-dimensional-benefit-risk-trade-off-matrix)
  - [Figure 12: Network Meta-Regression & Transitivity Bubble Plot](#figure-12-network-meta-regression--transitivity-bubble-plot)
- [Repository Architecture](#-repository-architecture)
- [Reproducibility & Execution Guide](#-reproducibility--execution-guide)
- [PRISMA-NMA Compliance Checklist](#-prisma-nma-compliance-checklist)
- [Citations & References](#-citations--references)

---

## 🌟 Executive Overview & Pedagogical Rationale

> [!NOTE]
> ### *"With this streamlined and reproducible evidence synthesis architecture, an investigator can build a complete evidence network, execute indirect comparisons, deconstruct multi-agent regimen synergy, and inspect a publication-grade League Table — grasping the core intuition of indirect evidence and transitivity before confronting the computational overhead and prior sensitivity of Bayesian MCMC modeling."*

When multiple interventions compete for the same clinical indication, head-to-head randomized controlled trials (RCTs) are rarely available for every pair of therapies. **Network Meta-Analysis (NMA)** resolves this therapeutic dilemma by synthesizing **direct evidence** (head-to-head trials) and **indirect evidence** (via common comparator arms) into a unified, coherent comparative framework.

This repository serves as a **portfolio-grade demonstration** of top-tier medical journal standards (*The Lancet*, *The New England Journal of Medicine*, *BMJ*, *JAMA Oncology*) for clinical trial synthesis, integrating:
1. **10 Production Statistical Analyses** (`scripts/analyses/`) covering estimation, ranking, inconsistency decomposition, influence cross-validation, component deconstruction, Monte Carlo simulation, benefit-risk modeling, and meta-regression.
2. **12 High-Resolution Figures (300 DPI)** (`scripts/designs/` & `outputs/figures/`) styled to perfection.
3. **Comprehensive HTML Report Dashboard** (`report/nma_comprehensive_report.html`) with embedded visual assets.

---

## 💡 Core Biostatistical & Methodological Competencies

This project systematically demonstrates ten advanced analytical capabilities essential for clinical biostatisticians, oncology researchers, and HTA evaluators:

1. **Systematic Survival Contrast Engineering ($\ln(\text{HR})$ and $\text{seTE}$):**
   Converting published survival endpoints (Hazard Ratios and $95\%$ Confidence Intervals) into symmetric Gaussian effect sizes:
   $$\text{TE} = \ln(\text{HR}), \quad \text{seTE} = \frac{\ln(\text{upper}) - \ln(\text{lower})}{2 \times 1.95996}$$

2. **Multi-Arm Trial Geometry & Correlation Handling:**
   Rigorous mathematical modeling of multi-arm trials enforcing **strict linear contrast additivity** ($\text{TE}_{BC} = \text{TE}_{AC} - \text{TE}_{AB}$) and positive-definite covariance matrices to satisfy `netmeta::chkmultiarm()`.

3. **Graph-Theoretical Network Inversion & Electrical Circuit Analogy:**
   Application of electrical circuit theory (Rücker 2012) where treatment nodes act as potentials $\mu_i$, trial comparisons act as conductive edges $w_{ij} = 1/\sigma_{ij}^2$, and network estimates are computed deterministically via the Moore-Penrose pseudoinverse ($L^+$) of the network Laplacian matrix.

4. **Treatment Ranking via Frequentist P-Scores & 10,000 Monte Carlo Draws:**
   Deterministic ranking of competing regimens measuring the certainty that one treatment is superior to another, paired with 10,000 multivariate normal draws yielding discrete rank probabilities ($P(\text{Rank} = r)$) and full cumulative rankograms (SUCRA).

5. **Leave-One-Out (LOO) Influence Cross-Validation:**
   Multi-threaded cross-validation sequentially omitting each of the 24 trials to quantify influence on the primary comparison (`IO + Chemo vs Chemo`), demonstrating remarkable stability (HR shift bounded within 0.033, 100% Rank 1 retention).

6. **Additive & Interactive Component Network Meta-Analysis (CNMA):**
   Deconstructing multi-agent oncology combinations into marginal active components (IO, CTLA4, TKI) and testing for pharmacologic synergy ($Q_{\text{diff}}$ test, $p = 0.0004$).

7. **Bi-dimensional Benefit-Risk Trade-Off Matrix:**
   Simultaneous dual NMA modeling mapping survival efficacy against severe Grade 3-5 toxicity in a 4-quadrant clinical decision matrix.

8. **Network Meta-Regression & Transitivity Diagnostics:**
   Formal evaluation of candidate effect modifiers (publication year, sample size $\ln(N)$, geographic setting) via `netmeta::netmetareg()`, confirming temporal stability ($\beta = 0.0003, p = 0.974$).

9. **Global & Local Inconsistency Diagnostics:**
   Orthogonal decomposition of Cochran's $Q$ ($Q = Q_{\text{het}} + Q_{\text{inc}}$), local node-splitting analysis (`netsplit`), and Net Heat matrix evaluation.

10. **Small-Study Effects & Publication Bias Evaluation:**
    Comparison-adjusted funnel plots ordered by established treatment hierarchy.

---

## 🎯 Clinical Research Scenario (Advanced NSCLC)

To ensure clinical relevance, this project synthesizes **First-Line Systemic Therapies for Advanced Non-Small Cell Lung Cancer (NSCLC)** without targetable driver mutations across **24 landmark Phase II/III Randomized Controlled Trials** covering **15,753 patients**:

| Regimen Code | Class & Mechanism | Representative Regimens | Role in Network |
|:---|:---|:---|:---|
| **`Chemo`** | Platinum-doublet Chemotherapy | Carboplatin/Cisplatin + Pemetrexed/Paclitaxel | **Standard Comparator (Anchor Reference)** |
| **`IO_Mono`** | Anti-PD-(L)1 Monotherapy | Pembrolizumab, Atezolizumab, Cemiplimab | Active Monotherapy |
| **`IO_Chemo`** | Immune Checkpoint Inhibitor + Chemotherapy | Pembrolizumab + Chemo, Tislelizumab + Chemo | Chemotherapy Combo |
| **`Dual_IO`** | Dual Checkpoint Blockade | Nivolumab + Ipilimumab (+/- limited chemo) | Chemotherapy-Free / Sparing |
| **`TKI`** | Targeted Tyrosine Kinase Inhibitor | Osimertinib, Gefitinib, Erlotinib | Targeted Monotherapy |
| **`TKI_Chemo`**| Targeted TKI + Platinum Chemotherapy | Osimertinib + Platinum/Pemetrexed | Targeted Combo |

---

## 🖼️ Publication Gallery (300 DPI Visual Exhibits)

| Figure | Description | High-Resolution Artifact |
|:---|:---|:---|
| **Figure 01** | Evidence Network Geometry (Weighted Topology) | [`outputs/figures/01_network_geometry.png`](outputs/figures/01_network_geometry.png) |
| **Figure 02** | Reference Comparison Forest Plot vs Chemotherapy | [`outputs/figures/02_forest_plot_random.png`](outputs/figures/02_forest_plot_random.png) |
| **Figure 03** | Frequentist P-Score Treatment Ranking Bar Chart | [`outputs/figures/03_pscore_ranking.png`](outputs/figures/03_pscore_ranking.png) |
| **Figure 04** | Node-Splitting Local Inconsistency Forest Plot | [`outputs/figures/04_netsplit_inconsistency.png`](outputs/figures/04_netsplit_inconsistency.png) |
| **Figure 05** | Net Heat Matrix Plot (Inconsistency & Hat Matrix) | [`outputs/figures/05_netheat_plot.png`](outputs/figures/05_netheat_plot.png) |
| **Figure 06** | Comparison-Adjusted Funnel Plot (Publication Bias) | [`outputs/figures/06_funnel_plot.png`](outputs/figures/06_funnel_plot.png) |
| **Figure 07** | Publication League Table Graphic Matrix (Random vs Direct) | [`outputs/figures/07_league_table_figure.png`](outputs/figures/07_league_table_figure.png) |
| **Figure 08** | Leave-One-Out (LOO) Sensitivity Forest Plot (24 Trials) | [`outputs/figures/08_leave_one_out_forest.png`](outputs/figures/08_leave_one_out_forest.png) |
| **Figure 09** | Component NMA (CNMA) Incremental Effects Forest | [`outputs/figures/09_component_effects.png`](outputs/figures/09_component_effects.png) |
| **Figure 10** | Probabilistic Hierarchy & Cumulative Rankograms | [`outputs/figures/10_rankograms.png`](outputs/figures/10_rankograms.png) |
| **Figure 11** | Bi-dimensional Benefit-Risk Trade-Off Matrix (OS vs Tox) | [`outputs/figures/11_benefit_risk_tradeoff.png`](outputs/figures/11_benefit_risk_tradeoff.png) |
| **Figure 12** | Network Meta-Regression Bubble & Transitivity Plot | [`outputs/figures/12_metaregression_bubble.png`](outputs/figures/12_metaregression_bubble.png) |

---

## 📁 Repository Architecture

```
nma-nsclc-evidence-synthesis/
├── data/
│   ├── nsclc_trial_contrasts.csv          # Clinical contrast dataset (34 contrasts, 24 RCTs)
│   └── nsclc_rob2_assessments.csv         # Risk of bias 2.0 evaluation data for web tools
├── scripts/
│   ├── analyses/
│   │   ├── 01_fit_nma_model.R             # Analysis 01: Model Estimation & Graph Laplacian Fit
│   │   ├── 02_treatment_rankings.R        # Analysis 02: P-Score Calculation & Hierarchy
│   │   ├── 03_league_table.R              # Analysis 03: Dual-Model League Table Generation
│   │   ├── 04_inconsistency_tests.R       # Analysis 04: Cochran's Q Global Decomposition
│   │   ├── 05_league_table_html.R         # Analysis 05: Formatted Interactive HTML League Table
│   │   ├── 06_leave_one_out_sensitivity.R # Analysis 06: Parallelized LOO Influence Cross-Validation
│   │   ├── 07_component_nma.R             # Analysis 07: Component NMA & Synergy Testing
│   │   ├── 08_rank_probabilities.R        # Analysis 08: 10,000 Monte Carlo Rank Probabilities & SUCRA
│   │   ├── 09_benefit_risk_tradeoff.R     # Analysis 09: Dual Efficacy vs Severe Toxicity NMA
│   │   └── 10_network_metaregression.R    # Analysis 10: Meta-Regression Across Year, Size & Region
│   ├── designs/
│   │   ├── fig01_network_geometry.R       # Design 01: Evidence Network Geometry (Topology)
│   │   ├── fig02_forest_plot.R            # Design 02: Reference Forest Plot vs Chemotherapy
│   │   ├── fig03_pscore_ranking.R         # Design 03: P-Score Ranking Bar Chart
│   │   ├── fig04_netsplit_inconsistency.R # Design 04: Node-Splitting Forest Plot
│   │   ├── fig05_netheat_plot.R           # Design 05: Net Heat Inconsistency Matrix
│   │   ├── fig06_funnel_plot.R            # Design 06: Comparison-Adjusted Funnel Plot
│   │   ├── fig07_league_table_matrix.R    # Design 07: Publication League Table Graphic Matrix
│   │   ├── fig08_leave_one_out_forest.R   # Design 08: Leave-One-Out Sensitivity Forest Plot
│   │   ├── fig09_component_effects.R      # Design 09: Component NMA Forest Plot
│   │   ├── fig10_rankograms.R             # Design 10: Multi-panel Rankograms & Cumulative Curves
│   │   ├── fig11_benefit_risk_tradeoff.R  # Design 11: 4-Quadrant Benefit-Risk Scatter Matrix
│   │   └── fig12_metaregression_bubble.R  # Design 12: Meta-Regression Bubble & Moderator Forest
│   └── run_all_pipeline.R                 # Master Orchestrator (10 Analyses + 12 Figures in ~35s)
├── outputs/
│   ├── figures/                           # 12 Publication-Grade 300 DPI PNG Figures
│   ├── tables/                            # Analysis CSV Tables + HTML League Table
│   └── models/                            # Serialized RDS Model Caches for Instant Downstream Builds
└── report/
    ├── nma_comprehensive_report.Rmd       # Comprehensive R Markdown source document
    └── nma_comprehensive_report.html      # Standalone 6.7 MB interactive HTML report
```

---

## 🚀 Reproducibility & Execution Guide

### Prerequisites
```r
install.packages(c("netmeta", "meta", "ggplot2", "readr", "knitr", "rmarkdown", "scales", "dplyr", "tidyr", "patchwork", "MASS", "parallel"))
```

### Execution Options

#### Option A: Run Full Master Pipeline (Analyses 1–10 + Figures 1–12)
```bash
Rscript scripts/run_all_pipeline.R
```
*Executes all 10 analyses and renders all 12 figures in ~35 seconds using smart RDS caching.*

#### Option B: Run Specific Individual Analysis or Design
```bash
# Example: Run Component NMA
Rscript scripts/analyses/07_component_nma.R
Rscript scripts/designs/fig09_component_effects.R

# Example: Run Benefit-Risk Trade-Off
Rscript scripts/analyses/09_benefit_risk_tradeoff.R
Rscript scripts/designs/fig11_benefit_risk_tradeoff.R
```

#### Option C: Render the Comprehensive HTML Report Dashboard
```r
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")
rmarkdown::render("report/nma_comprehensive_report.Rmd")
```

---

## ✅ PRISMA-NMA Compliance Checklist

| PRISMA-NMA Item | Guideline Description | Implementation in this Repository |
|:---|:---|:---|
| **Section 1: Title** | Identify report as a Network Meta-Analysis | Title explicitly states Frequentist Network Meta-Analysis |
| **Section 3: Rationale** | Explain need for indirect comparisons | Addressed in Theory section & League Table |
| **Section 6: Eligibility** | Define treatment nodes and trial criteria | 6 advanced NSCLC systemic regimens detailed |
| **Section 8: Geometry** | Present geometry of network graph | Figure 01: Weighted nodes (sample size) & edges (trials) |
| **Section 12: Synthesis** | Describe statistical methods for NMA | Graph-theoretical random/common effects via `netmeta` |
| **Section 14: Inconsistency** | Methods to assess consistency/transitivity | Global $Q$ test + Node-splitting (`netsplit`) + Net Heat (Figure 05) |
| **Section 16: Publication Bias** | Methods to evaluate small-study effects | Figure 06: Comparison-adjusted funnel plot |
| **Section 21: Study Results** | Provide summary data for each trial | Contrast dataset in `data/nsclc_trial_contrasts.csv` |
| **Section 22: Sensitivity** | Evaluate stability across trials | Figure 08: Leave-one-out cross-validation across 24 RCTs |
| **Section 23: Synthesis Results** | Present League Tables & treatment rankings | Figure 07 (League Table), Figure 03 (P-scores), Figure 10 (Rankograms) |
| **Section S1: Component NMA** | Deconstruct multi-agent combinations | Figure 09: Additive & interactive CNMA via `netcomb` |
| **Section S2: Benefit-Risk** | Multidimensional efficacy vs toxicity | Figure 11: 4-Quadrant survival vs Grade 3–5 severe toxicity matrix |
| **Section S3: Meta-Regression** | Screen candidate effect modifiers | Figure 12: Meta-regression across time, sample size, and region |

---

## 📚 Citations & References

1. **Rücker, G.** (2012). *Network meta-analysis, electrical networks and graph theory.* **Research Synthesis Methods**, 3(4), 312–324.
2. **Rücker, G., Petropoulou, M., & Schwarzer, G.** (2020). *Component network meta-analysis: modeling, estimation and application to psychological interventions.* **Biostatistics**, 21(4), 808–824.
3. **Salanti, G., et al.** (2011). *Evaluating the quality of evidence from a network meta-analysis.* **PLoS ONE**, 9(7), e99682.
4. **Jansen, J. P., & Naci, H.** (2013). *Conducting indirect-treatment-comparison and network-meta-analysis studies.* **Value in Health**, 14(4), 429–436.
5. **Hutton, B., et al.** (2015). *The PRISMA extension statement for reporting of systematic reviews incorporating network meta-analyses of health care interventions.* **Annals of Internal Medicine**, 162(11), 777–784.

---

<p align="center">
  <b>Engineered with scientific precision and methodological rigor.</b><br>
  <i>Designed for top-tier peer review and clinical decision support.</i>
</p>
