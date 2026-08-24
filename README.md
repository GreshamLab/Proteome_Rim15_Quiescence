# Proteome_Rim15_Quiescence

Analysis code for:

> Sun S. et al. "Rim15 orchestrates mitochondrial and cytoplasmic remodeling in
> quiescent *Saccharomyces cerevisiae* cells." PLOS Genetics (in revision).

SILAC-based proteomics and phosphoproteomics of wild-type and *rim15Δ* yeast
entering quiescence under carbon (C), nitrogen (N), and phosphorus (P)
starvation, sampled at 0, 6, 16, and 30 h, combined with physiological
characterization (growth, viability, bud index, DNA content), mitochondrial
imaging, and GEM nanoparticle single-particle tracking.

## Data availability

- Mass-spectrometry proteomics data: ProteomeXchange **PXD044239**.
- MaxQuant summary tables used as pipeline inputs are included in `data/`:
  `proteinGroups.csv` and `Phospho (STY)Sites.csv`.

## Repository layout

| Path | Contents |
|---|---|
| `data/` | MaxQuant input tables |
| `figure1/` | Self-contained code + source data for Figure 1 (growth, bud index, viability, SILAC label incorporation) |
| `tables/` | Intermediate and output tables produced by the scripts |
| `plots/` | Output figures |
| `output/wgcna/` | WGCNA network objects and module plots |
| `Dan/` | Independent normalization/QC check of the SILAC protein data |
| `rescue/drive-versions/` | Snapshots of script versions recovered from the lab's Google Drive working copy, kept pending reconciliation |

## Analysis pipeline

Scripts are numbered in run order and are written to be executed with the
repository root as the working directory (open `Proteome_Rim15_Quiescence.Rproj`
in RStudio and knit, or `rmarkdown::render()` from the root).

1. `01_data_preparation_qc.Rmd` — reads `data/proteinGroups.csv`; SILAC ratio
   normalization, QC, PCA, volcano plots; writes master protein tables to
   `tables/`.
2. `02_protein_anova_wgcna.Rmd` — two-way ANOVA of protein abundance;
   WGCNA co-expression modules per nutrient × genotype; KEGG module maps
   (outputs in `output/wgcna/`).
3. `03_ancova_analysis.Rmd` — ANCOVA on the tidy master protein table.
4. `04_protein_de_go_analysis.Rmd` — differential-expression summaries and
   GO/KEGG enrichment visualization.
5. `05_psites_analysis.Rmd` — phosphosite-level analysis from
   `data/Phospho (STY)Sites.csv`: normalization, differential phosphorylation,
   motif enrichment (rmotifx).
6. `06_pro_phospho_combined.Rmd` — joint proteome/phosphoproteome comparisons.
7. `07_single_protein_check.Rmd` — utility for inspecting individual proteins.

`figure1/Figure1.Rmd` is self-contained (run with `figure1/` as the working
directory).

Former script names (pre-2026 tidy): `Data preperation and QC.Rmd`,
`Proteom analysis - ANOVA.Rmd`, `ANCOVA analysis.Rmd`, `Proteom analysis.Rmd`,
`pSites analysis.Rmd`, `Pro-phospho combined.Rmd`, `Single protein check.Rmd`.

## R dependencies

Core: tidyverse (dplyr, tidyr, readr, ggplot2), reshape2, broom.
Statistics/omics: limma, qvalue, WGCNA, clusterProfiler, DOSE,
org.Sc.sgd.db, rmotifx, psych.
Visualization: pheatmap, gplots, ggpubr, ggrepel, ggfortify, factoextra,
cowplot, gridExtra, RColorBrewer, scales, UpSetR, plotly, GGally, ggformula.

Bioconductor packages (limma, qvalue, clusterProfiler, DOSE, org.Sc.sgd.db)
install via `BiocManager::install()`; rmotifx installs from GitHub
(`devtools::install_github("omarwagih/rmotifx")`).
