# Proteome_Rim15_Quiescence

Analysis code for:

> Sun S. et al. "Rim15 orchestrates mitochondrial and cytoplasmic remodeling in
> quiescent *Saccharomyces cerevisiae* cells." PLOS Genetics (accepted 2026).

SILAC-based proteomics and phosphoproteomics of wild-type and *rim15Δ* yeast
entering quiescence under carbon (C), nitrogen (N), and phosphorus (P)
starvation, sampled at 0, 6, 16, and 30 h, combined with physiological
characterization (growth, viability, bud index, DNA content), mitochondrial
imaging, and GEM nanoparticle single-particle tracking.

## Data availability

- Mass-spectrometry proteomics data: ProteomeXchange **PXD044239**.
- MaxQuant summary tables used as pipeline inputs are included in `data/`:
  `proteinGroups.csv` and `Phospho (STY)Sites.csv`.
- Reference datasets for cross-study comparisons are in
  `data/Publications for comparative analysis/` (see `data/Publications for comparative analysis/README.md`).
- Mitochondrial subcompartment annotation: `data/mitochondrial proteome/mmc4.xlsx`
  (Morgenstern et al. 2017, *Cell* supplementary table).
- Mitochondrial imaging and GEM tracking raw data: Mendeley Data
  (doi: [INSERT DOI — Figs 5–6]).
- **Source data:** Per-panel numerical data files are in `figureN/source_data/`
  and `sfig/source_data/`, provided per PLOS Genetics data deposition policy.

## Repository layout

| Path | Contents |
|---|---|
| `data/` | MaxQuant input tables and reference datasets |
| `figure1/` – `figure6/` | Self-contained Rmd + `source_data/` + `plots/` per figure |
| `sfig/` | Supplemental figure Rmds + `source_data/` + `plots/` |
| `analyses/` | Flow cytometry, cross-study comparisons, and imaging analyses |
| `tables/` | Intermediate and final output tables |
| `output/wgcna/` | WGCNA network objects and module plots |
| `Final_Supplemental_Tables/` | Assembled supplemental tables S1–S28 |

## How to run

### Prerequisites

Open `Proteome_Rim15_Quiescence.Rproj` in RStudio. All scripts use the **repo root**
as working directory.

### R package dependencies

Core: `tidyverse`, `reshape2`, `broom`.

Statistics/omics: `limma`, `qvalue`, `WGCNA`, `clusterProfiler`, `DOSE`,
`org.Sc.sgd.db`, `rmotifx`, `psych`.

Visualization: `pheatmap`, `gplots`, `ggpubr`, `ggrepel`, `ggfortify`,
`factoextra`, `cowplot`, `gridExtra`, `RColorBrewer`, `scales`, `UpSetR`,
`plotly`, `GGally`, `ggformula`, `ggridges`.

Bioconductor packages install via `BiocManager::install()`.
`rmotifx` installs from GitHub:
`devtools::install_github("omarwagih/rmotifx")`.

For Figure 4G: `install.packages("STRINGdb")`.

### Execution order

#### Step 1 — Main pipeline (run in order)

```r
rmarkdown::render("01_data_preparation_qc.Rmd")
rmarkdown::render("02_protein_anova_wgcna.Rmd")  # computationally intensive (~2 h)
source("extract_wgcna_tables.R")                 # exports WGCNA CSVs
rmarkdown::render("04_protein_de_go_analysis.Rmd")
rmarkdown::render("05_psites_analysis.Rmd")
rmarkdown::render("06_pro_phospho_combined.Rmd")
```

#### Step 2 — Additional analyses (run after Step 1; some have prerequisites)

```r
# Flow cytometry (requires interactive RStudio for CytoExploreR gating):
rmarkdown::render("analyses/flow_cytometry_gating.Rmd")  # prerequisite for SFig_7_8
rmarkdown::render("analyses/flow_ridge_plots.Rmd")

# Cross-study comparisons:
rmarkdown::render("analyses/phosphosite_crossstudy.Rmd") # prerequisite for SFig_10_11
rmarkdown::render("analyses/proteome_crossstudy.Rmd")
rmarkdown::render("analyses/mito_subcompartment.Rmd")    # requires mmc4.xlsx (see Data)
rmarkdown::render("analyses/combined_pca.Rmd")
```

#### Step 3 — Figure scripts (run after Steps 1–2)

```r
rmarkdown::render("figure1/Figure1.Rmd")
rmarkdown::render("figure2/Figure2.Rmd")
rmarkdown::render("figure3/Figure3.Rmd")
rmarkdown::render("figure4/Figure4.Rmd")
# figure5 and figure6: panels blocked pending collaborator data (see below)
rmarkdown::render("sfig/SFig_1_2.Rmd")
rmarkdown::render("sfig/SFig_3.Rmd")
rmarkdown::render("sfig/SFig_4.Rmd")
rmarkdown::render("sfig/SFig_5.Rmd")
rmarkdown::render("sfig/SFig_7_8.Rmd")   # requires flow_cytometry_gating.Rmd first
rmarkdown::render("sfig/SFig_9.Rmd")     # requires mmc4.xlsx (see Data)
rmarkdown::render("sfig/SFig_10_11.Rmd") # requires phosphosite_crossstudy.Rmd first
# SFig_6: blocked pending collaborator data (see below)
```

#### Step 4 — Supplemental tables

```r
source("make_supplemental_tables.R")
source("make_final_supplemental_tables.R")
```

### Figures 5 and 6, and S6 Fig (imaging data — BLOCKED)

Panels requiring per-cell segmentation and GEM trajectory data from Ying Xie
are marked **BLOCKED** in the respective Rmd files. Once the data is received:

- Mitochondrial volume (Fig 5B, S6 Fig): run `analyses/mito_volume.Rmd`, then
  uncomment the placeholder chunks in `figure5/Figure5.Rmd`.
- GEM diffusion (Fig 6B–D): run `analyses/gem_diffusion.Rmd`, then uncomment
  the placeholder chunks in `figure6/Figure6.Rmd`.

## Analysis pipeline (numbered scripts)

### Main pipeline

| Script | Inputs | Key outputs |
|---|---|---|
| `01_data_preparation_qc.Rmd` | `data/proteinGroups.csv` | `tables/master_tidy_protein.csv`, `tables/normed_ratio_perseus.csv`, QC/PCA plots |
| `02_protein_anova_wgcna.Rmd` | `tables/master_tidy_protein.csv` | `tables/anova_genotype.csv`, `tables/anova_nutrient.csv`, `tables/anova_genotype_nutrient.csv`, WGCNA RData in `output/wgcna/` |
| `04_protein_de_go_analysis.Rmd` | `tables/master_tidy_protein.csv`, ANOVA CSVs | `tables/DE_proteins_byGenotype.csv`, `tables/DE_proteins_byNutrient.csv`, GO/KEGG plots |
| `05_psites_analysis.Rmd` | `data/Phospho (STY)Sites.csv` | `tables/normed_pSite_ratio_perseus.csv`, `tables/DE_pSites_byGenotype.csv`, motif tables |
| `06_pro_phospho_combined.Rmd` | `tables/DE_proteins_byGenotype.csv`, `tables/DE_pSites_byGenotype.csv` | Combined proteome/phosphoproteome plots |
| `figure1/Figure1.Rmd` | `figure1/data/` | Figure 1 panels; run with `figure1/` as working directory |

### Helper scripts

| Script | Purpose | Inputs | Outputs |
|---|---|---|---|
| `extract_wgcna_tables.R` | Export WGCNA module assignments | `output/wgcna/*.RData` | `tables/wgcna_module_membership.csv`, `tables/wgcna_module_eigengenes.csv` |
| `make_supplemental_tables.R` | Assemble supplemental tables S5–S8 | `tables/` CSVs | `tables/supplemental/S5–S8_*.csv` |
| `make_final_supplemental_tables.R` | Assemble supplemental tables S9–S28 | multiple | `Final_Supplemental_Tables/` |

### Additional analyses (`analyses/`)

| Script | Purpose | Key outputs |
|---|---|---|
| `flow_cytometry_gating.Rmd` | Flow cytometry gating (interactive) | `analyses/output/dna_b2a_raw.csv`, `analyses/output/fsc_raw.csv` |
| `flow_ridge_plots.Rmd` | DNA content and cell size ridge plots | Input for SFig_7_8 |
| `phosphosite_crossstudy.Rmd` | Phosphosite cross-dataset overlap | `analyses/output/phosphosite_membership.csv`, `analyses/output/phosphoprotein_membership.csv` |
| `mito_subcompartment.Rmd` | Mitochondrial subcompartment trajectories | S9 Fig plots |
| `proteome_crossstudy.Rmd` | Proteome cross-dataset overlap | `analyses/output/proteome_coverage.csv` |
| `combined_pca.Rmd` | Combined proteome + phosphoproteome PCA | `analyses/output/combined_pca.pdf` |
| `mito_volume.Rmd` | Mito volume normalization (**BLOCKED**) | Fig 5B, S6 Fig |
| `gem_diffusion.Rmd` | GEM monomer-gated diffusion (**BLOCKED**) | Fig 6B–D |

## Figure–script mapping

| Figure  | Script                              | Panels (R-generated)                          | Notes |
|---------|-------------------------------------|-----------------------------------------------|-------|
| Fig 1   | `figure1/Figure1.Rmd`               | B (AA incorporation), C (growth), D (growth + bud index), E (1N DNA), F (viability) | Panel A = Illustrator diagram |
| Fig 2   | `figure2/Figure2.Rmd`               | A (quiescence protein heatmap), B (mito trajectories) | Sources: `04_protein_de_go_analysis.Rmd`, `01_data_preparation_qc.Rmd` |
| Fig 3   | `figure3/Figure3.Rmd`               | B (WPM eigengene heatmap), C (WPM × variable correlation), D (GO/KEGG), E (per-WPM distributions) | Panel A = Illustrator diagram |
| Fig 4   | `figure4/Figure4.Rmd`               | A, B, C, D, E, F, G                           | Sources: `05_psites_analysis.Rmd`; Panel G network layout in Cytoscape |
| Fig 5   | `figure5/Figure5.Rmd`               | B (mito volume; **BLOCKED**)                  | Panel A = microscopy images; Panel B awaiting data from Ying Xie |
| Fig 6   | `figure6/Figure6.Rmd`               | B, C, D (GEM; **BLOCKED**)                    | Panel A = microscopy images; B–D awaiting data from Ying Xie |
| S1 Fig  | `sfig/SFig_1_2.Rmd`                 | A (label incorporation), B (scan counts)      | |
| S2 Fig  | `sfig/SFig_1_2.Rmd`                 | A (replicate scatter), B (PCA)               | |
| S3 Fig  | `sfig/SFig_3.Rmd`                   | A (ANCOVA volcano), B (DE counts)             | |
| S4 Fig  | `sfig/SFig_4.Rmd`                   | A (phospho vs protein FC), B (target protein levels) | |
| S5 Fig  | `sfig/SFig_5.Rmd`                   | A (diagram), B (PEM heatmap), C (GO), D (dynamics) | Panel A = Illustrator |
| S6 Fig  | `sfig/SFig_6.Rmd` (**BLOCKED**)     | —                                             | Awaiting per-cell segmentation data from Ying Xie |
| S7 Fig  | `sfig/SFig_7_8.Rmd`                 | DNA ridge plots (WT + rim15Δ, C/N/P)         | Requires `analyses/flow_cytometry_gating.Rmd` run first |
| S8 Fig  | `sfig/SFig_7_8.Rmd`                 | FSC ridge plots (WT + rim15Δ, C/N/P)         | Requires `analyses/flow_cytometry_gating.Rmd` run first |
| S9 Fig  | `sfig/SFig_9.Rmd`                   | Mito subcompartment trajectories              | Requires `data/mitochondrial proteome/mmc4.xlsx` |
| S10 Fig | `sfig/SFig_10_11.Rmd`               | Phosphosite UpSet plot                        | Requires `analyses/phosphosite_crossstudy.Rmd` run first |
| S11 Fig | `sfig/SFig_10_11.Rmd`               | Phosphoprotein UpSet plot                     | Requires `analyses/phosphosite_crossstudy.Rmd` run first |

## Supplemental tables

| Table | File | Contents |
|---|---|---|
| S1 | `Final_Supplemental_Tables/Table S1.csv` | All protein SILAC quantification |
| S2 | `Final_Supplemental_Tables/Table S2.csv` | All phosphosite SILAC quantification |
| S3 | `tables/supplemental/S3_*.xlsx` | Proteins upregulated in C starvation required for quiescence |
| S4 | `tables/supplemental/S4_*.xlsx` | Proteins upregulated in P starvation required for quiescence |
| S5 | `tables/supplemental/S5_normalized_protein_ratios.csv` | Normalized SILAC ratios, 1,276 proteins × 74 samples |
| S6 | `tables/supplemental/S6_normalized_phosphosite_ratios.csv` | Normalized SILAC ratios, 5,056 phosphosites |
| S7 | `tables/supplemental/S7_ancova_results.csv` | ANCOVA results for all models |
| S8 | `tables/supplemental/S8_wgcna_module_membership.csv` | WGCNA module assignments |
| S9–S15 | `Final_Supplemental_Tables/` | Assembled by `make_final_supplemental_tables.R` |
| S16–S18 | `Final_Supplemental_Tables/` | WGCNA module tables |
| S19–S20 | `Final_Supplemental_Tables/` | Phosphosite DE results |
| S21–S22 | `Final_Supplemental_Tables/` | Cross-study phosphosite overlap |
| S23–S24 | `Final_Supplemental_Tables/` | Mitochondrial subcompartment tables |
| S25–S27 | `Final_Supplemental_Tables/` | Cross-study proteome overlap |
| S28 | `Final_Supplemental_Tables/` | Phosphosite membership |
