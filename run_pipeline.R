# run_pipeline.R
# Reproduces all main-text analyses and supplemental tables for Sun et al.
#
# Run from the repository root:
#   Rscript run_pipeline.R
#   — or —
#   source("run_pipeline.R")  # in RStudio with repo root as working directory
#
# Prerequisites:
#   - data/proteinGroups.csv and data/Phospho (STY)Sites.csv present
#   - output/wgcna/*.RData present (WGCNA is pre-computed; see Step 2 note)
#   - R packages listed in README.md installed
#
# Outputs: tables/, plots/, tables/supplemental/, figure1/output/
# Note: analyses/ scripts are run separately (flow_cytometry_gating.Rmd requires
#       an interactive RStudio session for CytoExploreR gating; other analyses/
#       scripts can be knitted non-interactively after gating is complete).

library(rmarkdown)

message("=== Step 1: Data preparation and QC ===")
render("01_data_preparation_qc.Rmd", quiet = TRUE)

# Step 2: WGCNA — computationally intensive (~hours per network).
# If output/wgcna/*.RData already exist (pre-computed), extract tables only.
# To rerun WGCNA from scratch, uncomment the render() line and comment out
# the source() line.
message("=== Step 2: WGCNA module table extraction ===")
# render("02_protein_anova_wgcna.Rmd", quiet = TRUE)  # rerun WGCNA (~hours)
source("extract_wgcna_tables.R")

# ANCOVA is run inside script 04.

message("=== Step 4: Differential abundance and GO/KEGG enrichment ===")
render("04_protein_de_go_analysis.Rmd", quiet = TRUE)

message("=== Step 5: Phosphosite analysis ===")
render("05_psites_analysis.Rmd", quiet = TRUE)

message("=== Step 6: Combined proteome + phosphoproteome ===")
render("06_pro_phospho_combined.Rmd", quiet = TRUE)

message("=== Figure 1 ===")
render("figure1/Figure1.Rmd",
       knit_root_dir = normalizePath("figure1"),
       quiet = TRUE)

message("=== Supplemental tables S5-S8 ===")
source("make_supplemental_tables.R")

message("\nPipeline complete.")
message("Main outputs: tables/, plots/, tables/supplemental/")
message("For additional analyses (flow cytometry, cross-study comparisons), see analyses/README.md")
