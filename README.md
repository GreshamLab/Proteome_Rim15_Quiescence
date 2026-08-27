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

Scripts are numbered in run order and should be executed with the repository
root as the working directory (open `Proteome_Rim15_Quiescence.Rproj` in
RStudio and knit, or `rmarkdown::render()` from the root).

### Main pipeline

| Script | Inputs | Key outputs |
|---|---|---|
| `01_data_preparation_qc.Rmd` | `data/proteinGroups.csv` | `tables/master_tidy_protein.csv`, `tables/normed_ratio_perseus.csv`, QC/PCA plots |
| `02_protein_anova_wgcna.Rmd` | `tables/master_tidy_protein.csv` | `tables/anova_genotype.csv`, `tables/anova_nutrient.csv`, `tables/anova_genotype_nutrient.csv`, WGCNA RData in `output/wgcna/` |
| `04_protein_de_go_analysis.Rmd` | `tables/master_tidy_protein.csv`, ANCOVA CSVs | `tables/DE_proteins_byGenotype.csv`, `tables/DE_proteins_byNutrient.csv`, GO/KEGG plots |
| `05_psites_analysis.Rmd` | `data/Phospho (STY)Sites.csv` | `tables/normed_pSite_ratio_perseus.csv`, `tables/DE_pSites_byGenotype.csv`, motif enrichment tables |
| `06_pro_phospho_combined.Rmd` | `tables/DE_proteins_byGenotype.csv`, `tables/DE_pSites_byGenotype.csv`, `tables/master_tidy_protein.csv`, ANCOVA CSVs | `plots/phospho_vs_protein_logFC.pdf` |
| `figure1/Figure1.Rmd` | `figure1/data/` | Figure 1 panels (growth curves, bud index, viability, label incorporation); run with `figure1/` as working directory |

### Helper scripts (run once, from repo root)

| Script | Purpose | Inputs | Outputs |
|---|---|---|---|
| `extract_wgcna_tables.R` | Export WGCNA module assignments from RData objects | `output/wgcna/*.RData` | `tables/wgcna_module_membership.csv`, `tables/wgcna_module_membership_wide.csv`, `tables/wgcna_module_eigengenes.csv` |
| `make_supplemental_tables.R` | Assemble submission-ready supplemental tables S5–S8 | `tables/master_tidy_protein.csv`, `tables/normed_pSite_ratio_perseus.csv`, ANCOVA CSVs, `tables/wgcna_module_membership_wide.csv` | `tables/supplemental/S5–S8_*.csv` |

### Third-revision analyses (`revision3/`)

Reviewer-required reanalyses for the third submission. All scripts use the
repository root as working directory. Outputs go to `revision3/output/`.

| Script | Reviewer point | Inputs | Key outputs |
|---|---|---|---|
| `R01_dna_content.Rmd` | R1.2 / R2.1 | FCS files in `data/Flow cytometry data/` (gitignored), `revision3/flow_sample_metadata.csv`, `revision3/flow_gating_template.csv` | `revision3/output/dna_b2a_raw.csv`, `revision3/output/fsc_raw.csv`, `dna_1n_percent_summary.csv`, `dna_1N_anova.csv` — **requires interactive RStudio session for CytoExploreR gating steps** |
| `R01b_ridge_plots.Rmd` | R1.2 / R2.1 | `revision3/output/dna_b2a_raw.csv`, `revision3/output/fsc_raw.csv` | `revision3/output/sfig_dna_ridges.pdf`, `revision3/output/sfig_fsc_ridges.pdf` |
| `R02_phosphosite_overlap.Rmd` | R1.3 / R2.5 | `tables/DE_pSites_byGenotype.csv`, published dataset CSVs | `revision3/output/phosphosite_overlap_upset.pdf`, `revision3/output/phosphosite_membership.csv` |
| `R03_mito_reconciliation.Rmd` | R2.7 | `tables/master_tidy_protein.csv`, Morgenstern 2017 reference list | `revision3/output/mito_coverage_table.csv`, `revision3/output/mito_by_class.pdf` |
| `R06_combined_pca.Rmd` | R2.6 | `tables/master_tidy_protein.csv` | `revision3/output/combined_pca.pdf` |

### Supplemental tables

Tables S1–S4 are in the journal submission package (`PLoS Genetics 3rd
Submission/Supplemental Information/` on Google Drive); S5–S8 are generated
by `make_supplemental_tables.R` and live in `tables/supplemental/`.

| Table | File | Contents |
|---|---|---|
| S1 | `Table S1.csv` | All protein SILAC quantification data in long format: intensity, normalized intensity, ratio, and normalized ratio for WT, rim15Δ, and spike-in across C/N/P starvation, timepoints 0/6/16/30 h, and all replicates (551,664 rows) |
| S2 | `Table S2.csv` | All phosphosite SILAC quantification data in long format: pSite ID, sequence window, genotype, nutrient, timepoint, replicate, and ratio value (106,200 rows) |
| S3 | `Table S3.xlsx` | The 11 proteins upregulated in carbon starvation (p.adj < 0.05) that are also genetically required for quiescence survival (p.adj < 0.05); columns: ORF, gene name, function |
| S4 | `Table S4.xlsx` | The 13 proteins upregulated in phosphorus starvation (p.adj < 0.05) that are also genetically required for quiescence survival; columns: ORF, gene name, function |
| S5 | `tables/supplemental/S5_normalized_protein_ratios.csv` | Normalized SILAC ratios for 1,276 proteins × 74 samples (WT + rim15Δ, C + P starvation, 0/6/16/30 h, 3 replicates), wide format |
| S6 | `tables/supplemental/S6_normalized_phosphosite_ratios.csv` | Normalized SILAC ratios for 5,056 phosphosites × 53 columns (full matrix), wide format |
| S7 | `tables/supplemental/S7_ancova_results.csv` | ANCOVA results for genotype, nutrient, and genotype × nutrient models combined |
| S8 | `tables/supplemental/S8_wgcna_module_membership.csv` | WGCNA module assignments for 927 proteins across 4 networks (C/P × WT/rim15Δ) |

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
