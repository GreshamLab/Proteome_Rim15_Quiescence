# Response letter text — R1.3 & R2.5 cross-dataset comparisons

Sun et al., PLOS Genetics PGENETICS-D-25-01150R1

Generated from analyses in `revision3/R02_phosphosite_overlap.Rmd` (phosphosite
overlap) and `revision3/R03_proteome_overlap.Rmd` (proteome coverage).

---

## R1.3 & R2.5 — Quantitative comparison to published phosphoproteome datasets

We thank both reviewers for this suggestion, which substantially strengthens the
contextualization of our results. We have now performed a comprehensive
cross-dataset comparison of our Rim15-regulated phosphosites against four
published Rim15/PP2A-Cdc55 phosphoproteomic datasets: Li et al. (2019), who
profiled the phosphoproteome of a *rim15Δ* strain under rich-medium conditions;
Dokladal et al. (2021), who identified high-confidence Rim15 substrates using a
thermal-shift proteomics approach; Baro et al. (2018), who catalogued
hyperphosphorylated peptides of putative PP2A-Cdc55 substrates; and Plank et al.
(2020), who profiled the phosphoproteome following acute auxin-mediated depletion
of the PP2A regulatory subunit Cdc55.

All datasets were harmonized to SGD systematic ORF names and absolute residue
positions in the canonical UniProt reference proteome (UP000002311), enabling
exact-site-level comparison. We validated the harmonization by confirming the
recovery of Igo1 S64 — the canonical, best-characterized Rim15 substrate site
(Talarek et al. 2010) — across multiple datasets. We present the overlap structure
as an UpSet plot (new Fig. SX) and provide a supplementary membership table
identifying every site and its dataset membership at two resolutions: exact site
(ORF + residue letter + position; Table S21) and protein level (Table S22).
Table S28 additionally provides side-by-side quantitative fold-change values from
each study for the overlapping sites.

We identified 30 phosphosites with significant Rim15-dependent regulation
(q ≤ 0.05) in carbon or phosphorus starvation, out of 4,941 phosphosites detected
in our experiment. Chi-square tests, using all detected sites as the universe,
showed that our Rim15-regulated sites were significantly enriched among those
reported by Dokladal et al. (2021; 3/30 sites overlap, 10.0% vs. 0.4% expected
from non-significant sites; χ² = 38.5, p = 5.6 × 10⁻¹⁰) and Baro et al. (2018;
9/30 sites, 30.0% vs. 10.9%; χ² = 9.3, p = 0.002). Overlap with Li et al. (2019;
1/30 sites, χ² = 0.3, p = 0.572) and Plank et al. (2020; 0/30 sites, p = 1.0)
was not statistically significant. Significant enrichment in two independent
datasets using distinct experimental strategies supports the validity of our
phosphoproteomic approach for identifying Rim15-dependent phosphorylation events.
The non-significant overlap with Li et al. and Plank et al. is not unexpected: Li
et al. used a different genetic background and profiled cells in rich medium rather
than starvation, and Plank et al.'s 10-minute Cdc55 depletion produced very few
robustly quantified exact sites (14 total), severely limiting power to detect
overlap with our starvation time-course data.

Of the 30 sites in our study, 10 (33%) appear in at least one published dataset
and 3 (10%) appear in two or more (consensus sites), including Igo1 S64 and sites
on the stress-responsive kinase targets Bcy1 and Sic1. The 20 sites unique to this
study represent candidate condition-specific Rim15 substrates under nutrient
starvation that were not observed in prior rich-medium profiling experiments. We
have rewritten the corresponding Discussion passage to reference this comparison
directly and to distinguish the consensus Rim15 targets from the condition-specific
sites identified here.

---

## Proteome coverage — supporting paragraph

*Addresses the representativeness of our proteome; can be incorporated into the
response to R2.7 (mitochondrial reconciliation) or the cover letter.*

To further contextualize the scope and representativeness of our proteome, we
compared proteins detected in our carbon-starvation WT time course against those
reported by Paulo et al. (2015), who used density-gradient sedimentation to
characterize the quiescent *S. cerevisiae* proteome by label-free mass
spectrometry. Of the 1,228 proteins we detected in WT cells during carbon
starvation, 1,101 (90%) were also detected by Paulo et al. (2015), who reported
2,437 proteins in total. The proteins detected only by Paulo et al. are consistent
with the greater dynamic range achievable by label-free quantification relative to
SILAC, which uses a common mid-log spike-in reference that constrains the effective
detection range. To assess directional consistency, we compared proteins regulated
≥2-fold at 30 h against the quiescence-enriched and log-enriched protein sets of
Paulo et al. (2015). Of proteins upregulated ≥2-fold in our study (n = 58), 16
(28%) were also enriched ≥2-fold in quiescent cells by Paulo et al. (χ² = 91.1,
p = 1.4 × 10⁻²¹); of proteins downregulated ≥2-fold (n = 23), 9 (39%) were
enriched in log-phase cells by Paulo et al. — that is, they decrease during
quiescence in both studies (χ² = 20.3, p = 6.7 × 10⁻⁶). This highly significant
directional agreement between independent experimental approaches demonstrates that
our SILAC-based proteome captures the major protein abundance transitions
associated with quiescence entry. A detailed quantitative comparison, including
Paulo et al. fold-change values where available, is provided in Tables S25–S27.
