# Additional analyses

Scripts for flow cytometry processing, cross-dataset comparisons, and
mitochondrial subcompartment analysis. All scripts run with the **repository
root** as working directory.

| Script | Purpose | Key outputs |
|---|---|---|
| `flow_cytometry_gating.Rmd` | Flow cytometry data import and CytoExploreR gating | `analyses/output/dna_b2a_raw.csv`, `analyses/output/fsc_raw.csv` — **requires interactive RStudio** |
| `flow_ridge_plots.Rmd` | Ridge plots of DNA content and cell size (run after gating) | Input for `sfig/SFig_7_8.Rmd` |
| `phosphosite_crossstudy.Rmd` | Cross-dataset phosphosite and phosphoprotein overlap | `analyses/output/phosphosite_membership.csv`, `analyses/output/phosphoprotein_membership.csv` — input for `sfig/SFig_10_11.Rmd` |
| `mito_subcompartment.Rmd` | Mitochondrial subcompartment trajectory analysis | `sfig/plots/SFig9_*.pdf` |
| `proteome_crossstudy.Rmd` | Cross-dataset proteome overlap | `analyses/output/proteome_coverage.csv`, `analyses/output/go_*_overlap.csv` |
| `combined_pca.Rmd` | Combined proteome + phosphoproteome PCA | `analyses/output/combined_pca.pdf` |
| `mito_volume.Rmd` | Mitochondrial volume normalization to cell volume (**BLOCKED**) | Fig 5B, S6 Fig — awaiting per-cell segmentation data from Ying Xie |
| `gem_diffusion.Rmd` | GEM nanoparticle monomer-gated diffusion analysis (**BLOCKED**) | Fig 6B–D — awaiting GEM trajectory data from Ying Xie |

Outputs are written to `analyses/output/`.
