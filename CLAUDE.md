# Proteome_Rim15_Quiescence

Analysis code for Sun et al., "Rim15 orchestrates mitochondrial and cytoplasmic
remodeling in quiescent Saccharomyces cerevisiae cells" (PLOS Genetics,
PGENETICS-D-25-01150R1, third revision in progress as of Aug 2026).

SILAC proteomics + phosphoproteomics of WT vs rim15Δ yeast entering quiescence
under carbon (C), nitrogen (N), and phosphorus (P) starvation, sampled at
0/6/16/30 h. Proteomics data: ProteomeXchange PXD044239.

## Working rules

- All analysis code is R / RMarkdown. Keep new analyses in R for consistency.
- Never run git operations in the Google Drive copy of this repo
  (`…/Sun et al.,/scripts/Proteome_Rim15_Quiescence`) — it is retired; Drive
  sync corrupted its index. This local clone (`~/Projects/Proteome_Rim15_Quiescence`)
  is the only working copy. Its unique content was rescued in commit 98ae7f1
  on branch `revision-3`.
- `rescue/drive-versions/` holds Drive copies of three Rmds that diverge from
  the GitHub versions (pending reconciliation): Data preperation and QC,
  Proteom analysis, pSites analysis. The Drive pSites version has post-2020
  fixes (VROOM_CONNECTION_SIZE, drop_na); the GitHub Data-prep version is the
  polished one. Reconcile before relying on either.

## External locations (Google Drive)

Base: `/Users/david/Library/CloudStorage/GoogleDrive-dg107@nyu.edu/My Drive/Gresham Lab_Papers/2026/Sun et al.,/`

- `PLoS Genetics 3rd Submission/` — decision letter PDF, response plan
  (`Sun_et_al_response_plan.docx` — the authoritative revision work plan with
  per-reviewer-point strategy, draft responses, and analysis recipes), final
  figures (PDF/TIFF), supplemental tables S1–S3.
- `data/` — MaxQuant Output, raw MS files, PRIDE submission package.
- Manuscript and point-by-point letter are Google Docs in the 3rd Submission
  folder; text revisions are delivered as a change-list (old → new text) that
  David pastes into Docs with Suggesting mode on. Do not edit the Docs directly.

## Script → figure map (to be completed during tidy)

- `Figure 1/Figure1.Rmd` — growth curves, bud index, viability (Fig 1).
- `Data preperation and QC.Rmd` — SILAC ratio prep, normalization, QC, PCA (Fig S2).
- `Proteom analysis.Rmd` / `Proteom analysis - ANOVA.Rmd` / `ANCOVA analysis.Rmd`
  — differential protein abundance (Figs 2–3).
- `pSites analysis.Rmd` — phosphosite analysis, motif enrichment (Fig 4).
- `Pro-phospho combined.Rmd` — combined proteome/phospho analyses.
- `Single protein check.Rmd` — per-protein inspection utility.
- WGCNA outputs: `*_networkConstruction-auto.RData`, `ME_*.pdf`, `*_WPM_dynamic*.pdf`.

## Third-revision analysis checklist (reviewer-required; recipes in response plan)

1. DNA-content reanalysis: gating strategy, 1N/2N histograms per genotype ×
   condition × time, two-way ANOVA on 1N-fraction trajectories. Needs FCS files (locate).
2. Cross-dataset phosphosite comparison (UpSet + table): Li 2019, Dokladal 2021,
   Baro 2018, Plank 2021; harmonize to ORF + residue.
3. Mitochondrial proteome reconciliation: coverage vs curated mito reference
   (Morgenstern 2017 / GO:0005739), WT vs rim15Δ by subcompartment/class.
4. Mito volume normalized to cell volume from cellpose masks. Needs masks + MitoGraph output (locate).
5. GEM monomer-intensity-gated Deff reanalysis. Needs trajectory data (locate).
6. Combined-genotype PCA (WT + rim15Δ), inspect WT carbon outlier replicate.
7. Numerical source data per graph + BioImage Archive deposition of imaging data.
