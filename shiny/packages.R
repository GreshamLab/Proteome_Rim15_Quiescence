# packages.R — install all Shiny app dependencies
# Run once: source("shiny/packages.R")

install.packages(c(
  "shiny",
  "bslib",
  "tidyverse",
  "plotly",
  "DT",
  "shinycssloaders"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "AnnotationDbi",
  "org.Sc.sgd.db",
  "GO.db"
))
