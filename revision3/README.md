# Third-revision analyses (PLOS Genetics, PGENETICS-D-25-01150R1)

Reviewer-required reanalyses for the third submission. Full strategy, draft
responses, and analysis recipes:
`…/Sun et al.,/PLoS Genetics 3rd Submission/Sun_et_al_response_plan.docx`.

All scripts run with the repository root as working directory.

| Script | Reviewer point | Status | Blocking data |
|---|---|---|---|
| `R01_dna_content.Rmd` | R1.2 / R2.1 | blocked | FCS files + gating info — **not in Drive or repo; ask Siyu** |
| `R02_phosphosite_overlap.Rmd` | R1.3 / R2.5 | ready to run | published supplements in `data/published/` (local, gitignored) |
| `R03_mito_reconciliation.Rmd` | R2.7 | ready to run | needs mito reference table (Morgenstern 2017 supp; download step in script) |
| `R04_mito_volume_norm.Rmd` | R2.8 | blocked | cellpose masks + MitoGraph output — **not in Drive; ask Siyu** |
| `R05_gem_monomer_gated.Rmd` | R2.9 | blocked | GEM trajectory + per-particle intensity data — **not in Drive; ask Siyu** |
| `R06_combined_pca.Rmd` | R2.6 | ready to run | none (uses `tables/normed_ratio_perseus.csv`) |
| `R07_source_data_deposition.md` | R2.12 | checklist | imaging raw data for BioImage Archive upload |

Outputs go to `revision3/output/` (figures + tables destined for new
supplementary figures/tables and the point-by-point response).
