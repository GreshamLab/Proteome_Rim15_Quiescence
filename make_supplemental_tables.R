# make_supplemental_tables.R
# Run from repo root:  Rscript make_supplemental_tables.R
#
# Produces four supplemental tables in tables/supplemental/:
#   S4_normalized_protein_ratios.csv
#   S5_normalized_phosphosite_ratios.csv
#   S6_ancova_results.csv
#   S7_wgcna_module_membership.csv

suppressPackageStartupMessages(library(tidyverse))

dir.create("tables/supplemental", showWarnings = FALSE)

# ── S4: Normalized protein ratios ────────────────────────────────────────────
# Source: master_tidy_protein.csv (long format, one row per protein × sample)
# Output: wide — one row per protein, columns = genotype_nutrient_time_rep
cat("Writing S4...\n")

protein_long <- read_csv("tables/master_tidy_protein.csv",
                         show_col_types = FALSE) %>%
  filter(genotype %in% c("WT", "rim15KO"),
         !is.na(normalized_ratio)) %>%
  mutate(sample = paste(genotype, nutrient, samplingTime, replicate, sep = "_")) %>%
  select(`Protein IDs`, `Gene names`, sample, normalized_ratio)

s4 <- protein_long %>%
  pivot_wider(names_from = sample, values_from = normalized_ratio) %>%
  arrange(`Gene names`)

write_csv(s4, "tables/supplemental/S4_normalized_protein_ratios.csv")
cat(sprintf("  S4: %d proteins × %d columns\n", nrow(s4), ncol(s4)))

# ── S5: Normalized phosphosite ratios ────────────────────────────────────────
# Source: normed_pSite_ratio_perseus.csv (already wide, one row per pSite)
# Drop: id, Fasta headers, WT/rim15KO ratio columns (keep WT and rim15KO only)
cat("Writing S5...\n")

psite_raw <- read_csv("tables/normed_pSite_ratio_perseus.csv",
                      show_col_types = FALSE)

# Keep metadata columns + WT and rim15KO sample columns (drop H/M ratio cols)
meta_cols  <- c("Proteins", "pSite", "Leading proteins", "Protein names", "Gene names")
# Sample columns: named {genotype}_{nutrient}_{time}_{rep}
# Drop the WT/rim15KO ratio channel (columns starting with "WT/rim15KO")
sample_cols <- names(psite_raw) %>%
  setdiff(c(meta_cols, "id", "Fasta headers")) %>%
  .[!startsWith(., "WT/rim15KO")]

s5 <- psite_raw %>%
  select(all_of(c(meta_cols, sample_cols))) %>%
  arrange(`Gene names`, pSite)

write_csv(s5, "tables/supplemental/S5_normalized_phosphosite_ratios.csv")
cat(sprintf("  S5: %d phosphosites × %d columns\n", nrow(s5), ncol(s5)))

# ── S6: ANCOVA results ────────────────────────────────────────────────────────
# Three models combined:
#   anova_genotype.csv          — lm(ratio ~ samplingTime * genotype) per nutrient
#   anova_nutrient.csv          — lm(ratio ~ samplingTime * nutrient) per genotype
#   anova_genotype_nutrient.csv — lm(ratio ~ samplingTime * genotype * nutrient)
# Filter to normed_ratio data; keep key terms only; add comparison column.
cat("Writing S6...\n")

# Key terms for each model
geno_terms    <- c("samplingTime", "genotypeWT", "samplingTime:genotypeWT")
nutri_terms   <- c("samplingTime", "nutrientP", "samplingTime:nutrientP")
threeway_terms <- c("samplingTime", "genotypeWT", "nutrientP",
                    "samplingTime:genotypeWT", "samplingTime:nutrientP",
                    "genotypeWT:nutrientP", "samplingTime:genotypeWT:nutrientP")

ancova_geno <- read_csv("tables/anova_genotype.csv", show_col_types = FALSE) %>%
  filter(data == "normed_ratio", term %in% geno_terms) %>%
  mutate(comparison = "genotype (per nutrient)",
         condition  = nutrient) %>%
  select(`Protein IDs`, `Gene names`, comparison, condition,
         term, estimate, std.error, statistic, p.value, qvalue)

ancova_nutri <- read_csv("tables/anova_nutrient.csv", show_col_types = FALSE) %>%
  filter(data == "normed_ratio", term %in% nutri_terms) %>%
  mutate(comparison = "nutrient (per genotype)",
         condition  = genotype) %>%
  select(`Protein IDs`, `Gene names`, comparison, condition,
         term, estimate, std.error, statistic, p.value, qvalue)

ancova_3way <- read_csv("tables/anova_genotype_nutrient.csv",
                        show_col_types = FALSE) %>%
  filter(data == "normed_ratio", term %in% threeway_terms) %>%
  mutate(comparison = "genotype × nutrient (combined)",
         condition  = NA_character_) %>%
  select(`Protein IDs`, `Gene names`, comparison, condition,
         term, df, statistic, p.value, qvalue)

# Bind — note anova_genotype_nutrient has df/sumsq not estimate/std.error
s6 <- bind_rows(ancova_geno, ancova_nutri, ancova_3way) %>%
  arrange(comparison, `Gene names`, term)

write_csv(s6, "tables/supplemental/S6_ancova_results.csv")
cat(sprintf("  S6: %d rows\n", nrow(s6)))

# ── S7: WGCNA module membership ───────────────────────────────────────────────
cat("Writing S7...\n")

s7 <- read_csv("tables/wgcna_module_membership_wide.csv",
               show_col_types = FALSE) %>%
  rename(`Protein IDs` = protein_id, `Gene names` = gene_name) %>%
  arrange(`Gene names`)

write_csv(s7, "tables/supplemental/S7_wgcna_module_membership.csv")
cat(sprintf("  S7: %d proteins\n", nrow(s7)))

cat("\nDone. Files in tables/supplemental/:\n")
cat("  S4_normalized_protein_ratios.csv\n")
cat("  S5_normalized_phosphosite_ratios.csv\n")
cat("  S6_ancova_results.csv\n")
cat("  S7_wgcna_module_membership.csv\n")
