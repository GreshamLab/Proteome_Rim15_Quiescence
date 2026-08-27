# Sun et al. — Revision 3 Work Plan

Paper: "Rim15 orchestrates mitochondrial and cytoplasmic remodeling in quiescent
Saccharomyces cerevisiae cells" (PLOS Genetics, PGENETICS-D-25-01150R1)


---

## Reproducibility & Documentation

- [x] Audit existing repo for scripts that already perform these analyses (protein ratio
      normalization, phosphosite normalization, ANCOVA, WGCNA, flow cytometry processing) —
      inventory what exists before writing anything new
- [x] Identify and label ambiguous/undocumented scripts — purpose/inputs/outputs blocks
      added to 01–06; scripts 03 and 07 marked RETIRED
- [x] Standardize output file naming — supplemental tables named S5–S8_*.csv;
      revision3 outputs named descriptively
- [x] Refactor pipeline so it reproducibly regenerates all main-text analyses and
      supplemental tables/figures — `run_pipeline.R` created as single entry point
- [x] Write/update README.md with script → input → output map

## New Supplemental Tables

- [x] S5 — Normalized protein SILAC ratios (`tables/supplemental/S5_normalized_protein_ratios.csv`): 1,276 proteins × 74 sample columns (WT + rim15Δ, C + P, T0/6/16/30, 3 reps)
- [x] S6 — Normalized phosphosite ratios (`tables/supplemental/S6_normalized_phosphosite_ratios.csv`): 5,056 pSites × 53 columns (full matrix, not filtered by significance)
- [x] S7 — ANCOVA results (`tables/supplemental/S7_ancova_results.csv`): genotype, nutrient, and 3-way models combined; filtered to normed_ratio data
- [x] S8 — WGCNA module membership (`tables/supplemental/S8_wgcna_module_membership.csv`): 927 proteins × 4 networks (C/P × WT/rim15Δ)

## New Supplemental Figures — Flow Cytometry

- [x] Ridge plots: DNA content over time, split by genotype (2) x condition (3) — `revision3/output/sfig_dna_ridges.pdf`
- [x] Cell size (FSC) over time, split by genotype (2) x condition (3) — `revision3/output/sfig_fsc_ridges.pdf`
- [x] Check for existing flow cytometry processing code before writing new

---


## Reviewer-required analyses

| ID | Script | Reviewer point | Status | Notes |
|---|---|---|---|---|
| R01 | `revision3/R01_dna_content.Rmd` | R1.2 / R2.1 | done | |
| R02 | `revision3/R02_phosphosite_overlap.Rmd` | R1.3 / R2.5 | done | |
| R03 | `revision3/R03_mito_reconciliation.Rmd` | R2.7 | done | |
| R06 | `revision3/R06_combined_pca.Rmd` | R2.6 | done | |
## Other tasks

- [ ] Shiny app deployment (NYU Research Computing)
- [x] rescue/drive-versions/ reconciliation — 01 and 04 repo versions are polished; 05 updated with drop_na(Proteins) and pSite[is.na(pSite)] <- 0 from Drive version

---

## Notes

