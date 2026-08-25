# Third-revision analyses (PLOS Genetics, PGENETICS-D-25-01150R1)

Reviewer-required reanalyses for the third submission. Full strategy, draft
responses, and analysis recipes:
`…/Sun et al.,/PLoS Genetics 3rd Submission/Sun_et_al_response_plan.docx`.

All scripts run with the repository root as working directory.

| Script | Reviewer point | Status | Blocking data |
|---|---|---|---|
| `R01_dna_content.Rmd` | R1.2 / R2.1 | ready to run | FCS files in `data/Flow cytometry data/` (gitignored); gating + metadata in `revision3/flow_*.csv` |
| `R02_phosphosite_overlap.Rmd` | R1.3 / R2.5 | **done** | — |
| `R03_mito_reconciliation.Rmd` | R2.7 | **done** | — |
| `R04_mito_volume_norm.Rmd` | R2.8 | blocked | cellpose masks + MitoGraph output — **not in Drive; ask Siyu** |
| `R05_gem_monomer_gated.Rmd` | R2.9 | blocked | GEM trajectory + per-particle intensity data — **not in Drive; ask Siyu** |
| `R06_combined_pca.Rmd` | R2.6 | **done** | — |
| `R07_source_data_deposition.md` | R2.12 | checklist | imaging raw data for BioImage Archive upload |

Outputs go to `revision3/output/` (figures + tables destined for new
supplementary figures/tables and the point-by-point response).
