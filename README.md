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
- Mitochondrial imaging and GEM tracking raw data: Mendeley Data
  (doi: [INSERT DOI — Figs 5–6]).
- **Source data:** Per-panel numerical data files are in `figureN/source_data/`
  and `sfig/source_data/`, provided per PLOS Genetics data deposition policy.

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
| S7 Fig  | `sfig/SFig_7_8.Rmd`                 | DNA ridge plots (WT + rim15Δ, C/N/P)         | Requires `revision3/R01_dna_content.Rmd` run first |
| S8 Fig  | `sfig/SFig_7_8.Rmd`                 | FSC ridge plots (WT + rim15Δ, C/N/P)         | Requires `revision3/R01_dna_content.Rmd` run first |
| S9 Fig  | `sfig/SFig_9.Rmd`                   | Mito subcompartment trajectories              | Requires `data/mitochondrial proteome/mmc4.xlsx` (local) |
| S10 Fig | `sfig/SFig_10_11.Rmd`               | Phosphosite UpSet plot                        | Requires `revision3/R02_phosphosite_overlap.Rmd` run first |
| S11 Fig | `sfig/SFig_10_11.Rmd`               | Phosphoprotein UpSet plot                     | Requires `revision3/R02_phosphosite_overlap.Rmd` run first |

## Repository layout

| Path | Contents |
|---|---|
| `data/` | MaxQuant input tables (`proteinGroups.csv`, `Phospho (STY)Sites.csv`) |
| `figure1/` | Self-contained code + source data for Figure 1 |
| `figure2/` | Figure 2 Rmd + `source_data/` + `plots/` |
| `figure3/` | Figure 3 Rmd + `source_data/` + `plots/` |
| `figure4/` | Figure 4 Rmd + `source_data/` + `plots/` |
| `figure5/` | Figure 5 Rmd stub (Panel B blocked) + `source_data/` + `plots/` |
| `figure6/` | Figure 6 Rmd stub (Panels B–D blocked) + `source_data/` + `plots/` |
| `sfig/` | Supplemental figure Rmds + `source_data/` + `plots/` |
| `tables/` | Intermediate and final output tables |
| `plots/` | Legacy output directory (numbered pipeline) |
| `output/wgcna/` | WGCNA network objects and module plots |
| `revision3/` | Third-revision analysis scripts (R01–R07) |
| `Dan/` | Independent QC check of SILAC protein data |

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

1. **Main pipeline** (must run in order; produces `tables/`):
   ```r
   rmarkdown::render("01_data_preparation_qc.Rmd")
   rmarkdown::render("02_protein_anova_wgcna.Rmd")  # computationally intensive
   source("extract_wgcna_tables.R")                 # exports WGCNA CSVs
   rmarkdown::render("04_protein_de_go_analysis.Rmd")
   rmarkdown::render("05_psites_analysis.Rmd")
   rmarkdown::render("06_pro_phospho_combined.Rmd")
   ```

2. **Third-revision analyses** (independent, can run after step 1):
   ```r
   # R01_dna_content.Rmd requires interactive RStudio for CytoExploreR gating
   rmarkdown::render("revision3/R01b_ridge_plots.Rmd")
   rmarkdown::render("revision3/R02_phosphosite_overlap.Rmd")
   rmarkdown::render("revision3/R03_mito_reconciliation.Rmd")
   ```

3. **Figure Rmds** (run after step 1; each is self-contained and reads from `tables/`):
   ```r
   rmarkdown::render("figure2/Figure2.Rmd")
   rmarkdown::render("figure3/Figure3.Rmd")
   rmarkdown::render("figure4/Figure4.Rmd")
   # figure5 and figure6 are stubs pending collaborator data
   rmarkdown::render("sfig/SFig_1_2.Rmd")
   rmarkdown::render("sfig/SFig_3.Rmd")
   rmarkdown::render("sfig/SFig_4.Rmd")
   rmarkdown::render("sfig/SFig_5.Rmd")
   rmarkdown::render("sfig/SFig_7_8.Rmd")  # requires revision3/R01b first
   rmarkdown::render("sfig/SFig_9.Rmd")    # requires local mmc4.xlsx
   rmarkdown::render("sfig/SFig_10_11.Rmd") # requires revision3/R02 first
   ```

### Figures 5 and 6 (imaging data)

Panels requiring collaborator data (Ying Xie) are marked **BLOCKED** in
`figure5/Figure5.Rmd` and `figure6/Figure6.Rmd`. Once the per-cell segmentation
masks and GEM trajectory data are received:
- Mitochondrial volume (Fig 5B, S6 Fig): run `revision3/R04_mito_volume_norm.Rmd`
- GEM diffusion (Fig 6B–D): run `revision3/R05_gem_monomer_gated.Rmd`
- Then uncomment the placeholder chunks in the respective figure Rmds.

## Analysis pipeline (numbered scripts)

### Main pipeline

| Script | Inputs | Key outputs |
|---|---|---|
| `01_data_preparation_qc.Rmd` | `data/proteinGroups.csv` | `tables/master_tidy_protein.csv`, `tables/normed_ratio_perseus.csv`, QC/PCA plots |
| `02_protein_anova_wgcna.Rmd` | `tables/master_tidy_protein.csv` | `tables/anova_genotype.csv`, `tables/anova_nutrient.csv`, `tables/anova_genotype_nutrient.csv`, WGCNA RData in `output/wgcna/` |
| `04_protein_de_go_analysis.Rmd` | `tables/master_tidy_protein.csv`, ANCOVA CSVs | `tables/DE_proteins_byGenotype.csv`, `tables/DE_proteins_byNutrient.csv`, GO/KEGG plots |
| `05_psites_analysis.Rmd` | `data/Phospho (STY)Sites.csv` | `tables/normed_pSite_ratio_perseus.csv`, `tables/DE_pSites_byGenotype.csv`, motif tables |
| `06_pro_phospho_combined.Rmd` | `tables/DE_proteins_byGenotype.csv`, `tables/DE_pSites_byGenotype.csv`, `tables/master_tidy_protein.csv`, ANCOVA CSVs | `plots/phospho_vs_protein_logFC.pdf` |
| `figure1/Figure1.Rmd` | `figure1/data/` | Figure 1 panels; run with `figure1/` as working directory |

### Helper scripts

| Script | Purpose | Inputs | Outputs |
|---|---|---|---|
| `extract_wgcna_tables.R` | Export WGCNA module assignments | `output/wgcna/*.RData` | `tables/wgcna_module_membership.csv`, `tables/wgcna_module_eigengenes.csv` |
| `make_supplemental_tables.R` | Assemble supplemental tables S5–S8 | `tables/` CSVs | `tables/supplemental/S5–S8_*.csv` |

### Third-revision analyses (`revision3/`)

| Script | Reviewer point | Key outputs |
|---|---|---|
| `R01_dna_content.Rmd` | R1.2/R2.1 | `revision3/output/dna_b2a_raw.csv`, `revision3/output/fsc_raw.csv` — **requires interactive RStudio** |
| `R01b_ridge_plots.Rmd` | R1.2/R2.1 | S7+S8 Fig ridge plots |
| `R02_phosphosite_overlap.Rmd` | R1.3/R2.5 | S10+S11 Fig UpSet plots + membership tables |
| `R03_mito_reconciliation.Rmd` | R2.7 | S9 Fig mito subcompartment trajectories |
| `R04_mito_volume_norm.Rmd` | R2.8 | Fig 5B, S6 Fig — **BLOCKED: data from Ying Xie** |
| `R05_gem_monomer_gated.Rmd` | R2.9 | Fig 6B–D — **BLOCKED: data from Ying Xie** |

## Supplemental tables

| Table | File | Contents |
|---|---|---|
| S1 | `Table S1.csv` | All protein SILAC quantification (551,664 rows) |
| S2 | `Table S2.csv` | All phosphosite SILAC quantification (106,200 rows) |
| S3 | `Table S3.xlsx` | 11 proteins upregulated in C starvation that are required for quiescence |
| S4 | `Table S4.xlsx` | 13 proteins upregulated in P starvation that are required for quiescence |
| S5 | `tables/supplemental/S5_normalized_protein_ratios.csv` | Normalized SILAC ratios, 1,276 proteins × 74 samples |
| S6 | `tables/supplemental/S6_normalized_phosphosite_ratios.csv` | Normalized SILAC ratios, 5,056 phosphosites |
| S7 | `tables/supplemental/S7_ancova_results.csv` | ANCOVA results for all models |
| S8 | `tables/supplemental/S8_wgcna_module_membership.csv` | WGCNA module assignments |

Former script names (pre-2026 tidy): `Data preperation and QC.Rmd`,
`Proteom analysis - ANOVA.Rmd`, `ANCOVA analysis.Rmd`, `Proteom analysis.Rmd`,
`pSites analysis.Rmd`, `Pro-phospho combined.Rmd`, `Single protein check.Rmd`.
