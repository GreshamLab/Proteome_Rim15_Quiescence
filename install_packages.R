# install_packages.R
# Run this script once before executing any pipeline scripts.
# Sets up all R package dependencies for Proteome_Rim15_Quiescence.
#
# Usage (from RStudio console or terminal):
#   source("install_packages.R")

# ── CRAN packages ─────────────────────────────────────────────────────────────
cran_pkgs <- c(
  "tidyverse",      # dplyr, tidyr, ggplot2, readr, etc.
  "reshape2",
  "broom",
  "psych",          # geometric mean
  "data.table",     # fast CSV reading in analyses/phosphosite_crossstudy.Rmd
  "readxl",
  "GGally",
  "ggpubr",
  "ggrepel",
  "ggfortify",
  "ggridges",
  "ggformula",
  "gplots",
  "pheatmap",
  "gridExtra",
  "cowplot",
  "RColorBrewer",
  "scales",
  "plotly",
  "UpSetR",
  "WGCNA",          # weighted gene co-expression network analysis
  "mosaicData",
  "SuperExactTest",
  "devtools",       # loaded by 01_data_preparation_qc.Rmd
  "remotes",        # for GitHub installs below
  "rmarkdown",
  "STRINGdb"        # Figure 4G protein interaction network
)

to_install <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) {
  message("Installing CRAN packages: ", paste(to_install, collapse = ", "))
  install.packages(to_install)
} else {
  message("All CRAN packages already installed.")
}

# ── Bioconductor packages ─────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_pkgs <- c(
  "limma",
  "qvalue",
  "clusterProfiler",
  "DOSE",
  "org.Sc.sgd.db",
  "GO.db",
  "AnnotationDbi",
  "Biostrings",
  # Flow cytometry (only needed for analyses/flow_cytometry_gating.Rmd)
  "cytolib",
  "flowCore",
  "flowWorkspace",
  "openCyto"
)

bioc_to_install <- bioc_pkgs[!bioc_pkgs %in% rownames(installed.packages())]
if (length(bioc_to_install) > 0) {
  message("Installing Bioconductor packages: ", paste(bioc_to_install, collapse = ", "))
  BiocManager::install(bioc_to_install, ask = FALSE)
} else {
  message("All Bioconductor packages already installed.")
}

# ── GitHub packages ───────────────────────────────────────────────────────────
# rmotifx — phosphorylation motif enrichment (used in 05_psites_analysis.Rmd)
if (!requireNamespace("rmotifx", quietly = TRUE)) {
  message("Installing rmotifx from GitHub (omarwagih/rmotifx)...")
  remotes::install_github("omarwagih/rmotifx")
} else {
  message("rmotifx already installed.")
}

# CytoExploreR — interactive flow cytometry gating
# (only needed for analyses/flow_cytometry_gating.Rmd)
if (!requireNamespace("CytoExploreR", quietly = TRUE)) {
  message("Installing CytoExploreR from GitHub (requires superheat first)...")
  if (!requireNamespace("superheat", quietly = TRUE))
    remotes::install_github("rlbarter/superheat")
  remotes::install_github("DillonHammill/CytoExploreR")
} else {
  message("CytoExploreR already installed.")
}

message("\nAll packages installed. You can now run the pipeline.")
