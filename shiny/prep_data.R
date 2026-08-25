# prep_data.R — run once from the repo root before launching the Shiny app:
#   Rscript shiny/prep_data.R
#
# Reads tables/ CSVs, trims to the columns the app needs, and saves compact
# .rds files to shiny/data/.  This step is done offline so the app never reads
# the raw CSVs at startup.
#
# Note: AnnotationDbi masks both dplyr::select and dplyr::rename, so every call
# to those functions uses an explicit dplyr:: prefix throughout this script.
# Columns with spaces/special chars in their names (e.g. "Gene names") are
# renamed via base R *before* entering the dplyr pipe, because backtick-quoted
# names in rename() / select() fail when those functions are masked.

suppressPackageStartupMessages({
  library(tidyverse)
  library(AnnotationDbi)   # loaded for GO lookup; masks dplyr::select/rename
  library(org.Sc.sgd.db)
  library(GO.db)           # needed for GO term descriptions (TERM column)
})

dir.create("shiny/data", recursive = TRUE, showWarnings = FALSE)

# ── 1. Protein abundance (long format) ─────────────────────────────────────
cat("Reading protein abundance...\n")
tmp <- read_csv("tables/master_tidy_protein.csv", show_col_types = FALSE)
# Rename space-containing column names via base R (NSE-safe)
names(tmp)[names(tmp) == "Gene names"]  <- "gene_raw"
names(tmp)[names(tmp) == "Protein IDs"] <- "protein_id"
# Now use dplyr::select() to pick and rename the remaining columns
# (avoids the duplicate "ratio" column that base-R renaming would cause)
protein_long <- tmp %>%
  dplyr::select(gene_raw, protein_id, genotype, nutrient,
                time       = samplingTime,
                replicate,
                ratio      = normalized_ratio) %>%
  mutate(
    gene     = str_extract(gene_raw, "^[^;]+"),   # first gene in protein group
    genotype = recode(genotype, "rim15KO" = "rim15δ"),
    time_h   = as.integer(str_extract(time, "\\d+"))
  ) %>%
  dplyr::select(gene, protein_id, genotype, nutrient, time_h, replicate, ratio) %>%
  filter(!is.na(ratio), ratio > 0, !is.na(gene))
rm(tmp)

saveRDS(protein_long, "shiny/data/protein_long.rds")
cat(sprintf("  protein_long.rds: %d rows, %d genes\n",
            nrow(protein_long), n_distinct(protein_long$gene)))

# ── 2. Protein differential abundance ──────────────────────────────────────
cat("Reading protein DE...\n")
tmp <- read_csv("tables/DE_proteins_byGenotype.csv", show_col_types = FALSE)
names(tmp)[names(tmp) == "Gene names"]  <- "gene_raw"
names(tmp)[names(tmp) == "Protein IDs"] <- "protein_id"

protein_de <- tmp %>%
  dplyr::select(gene_raw, protein_id, nutrient, time, pvalue, qvalue, log2FC) %>%
  mutate(gene = str_extract(gene_raw, "^[^;]+")) %>%
  dplyr::select(gene, protein_id, nutrient, time, pvalue, qvalue, log2FC) %>%
  filter(!is.na(log2FC), is.finite(log2FC), !is.na(gene))
rm(tmp)

saveRDS(protein_de, "shiny/data/protein_de.rds")
cat(sprintf("  protein_de.rds: %d rows\n", nrow(protein_de)))

# ── 3. Phosphosite abundance (pivot wide → long) ────────────────────────────
cat("Reading phosphosite abundance...\n")
psite_wide <- read_csv("tables/normed_pSite_ratio_perseus.csv",
                       show_col_types = FALSE)

meta_cols   <- c("Proteins", "id", "pSite", "Leading proteins",
                 "Protein names", "Gene names", "Fasta headers")
sample_cols <- setdiff(names(psite_wide), meta_cols)
# Keep only per-genotype channels; drop "WT/rim15KO" SILAC ratio columns
sample_cols <- sample_cols[!grepl("^WT/", sample_cols)]

tmp <- psite_wide %>%
  dplyr::select(all_of(c("pSite", "Gene names", sample_cols)))
names(tmp)[names(tmp) == "pSite"]      <- "psite"
names(tmp)[names(tmp) == "Gene names"] <- "gene_raw"

psite_long <- tmp %>%
  mutate(gene = str_extract(gene_raw, "^[^;]+")) %>%
  dplyr::select(-gene_raw) %>%
  pivot_longer(all_of(sample_cols), names_to = "sample", values_to = "ratio") %>%
  filter(!is.na(ratio), ratio > 0) %>%
  # sample column names: genotype_nutrient_timepoint_replicate
  separate(sample,
           into  = c("genotype", "nutrient", "time", "replicate"),
           sep   = "_",
           extra = "drop") %>%
  mutate(
    genotype = recode(genotype, "rim15KO" = "rim15δ"),
    time_h   = as.integer(str_extract(time, "\\d+"))
  ) %>%
  dplyr::select(psite, gene, genotype, nutrient, time_h, replicate, ratio) %>%
  filter(!is.na(gene))
rm(tmp, psite_wide)

saveRDS(psite_long, "shiny/data/psite_long.rds")
cat(sprintf("  psite_long.rds: %d rows, %d psites\n",
            nrow(psite_long), n_distinct(psite_long$psite)))

# ── 4. Phosphosite differential abundance ──────────────────────────────────
cat("Reading phosphosite DE...\n")
tmp <- read_csv("tables/DE_pSites_byGenotype.csv", show_col_types = FALSE)
names(tmp)[names(tmp) == "Gene names"] <- "gene_raw"
names(tmp)[names(tmp) == "pSite"]      <- "psite"

psite_de <- tmp %>%
  dplyr::select(gene_raw, psite, nutrient, time, pvalue, qvalue, log2FC) %>%
  mutate(gene = str_extract(gene_raw, "^[^;]+")) %>%
  dplyr::select(gene, psite, nutrient, time, pvalue, qvalue, log2FC) %>%
  filter(!is.na(log2FC), is.finite(log2FC), !is.na(gene))
rm(tmp)

saveRDS(psite_de, "shiny/data/psite_de.rds")
cat(sprintf("  psite_de.rds: %d rows\n", nrow(psite_de)))

# ── 5. Gene universe (for search + GO term expansion) ───────────────────────
cat("Building gene universe...\n")
genes_detected <- sort(unique(na.omit(protein_long$gene)))

# Map gene names → ORF via GENENAME keytype for GO term lookup
gene_orf <- tryCatch(
  AnnotationDbi::select(org.Sc.sgd.db,
                        keys    = genes_detected,
                        keytype = "GENENAME",
                        columns = "ORF") %>%
    filter(!is.na(ORF)) %>%
    distinct(GENENAME, ORF),
  error = function(e) {
    message("  GENENAME keytype failed, trying COMMON")
    AnnotationDbi::select(org.Sc.sgd.db,
                          keys    = genes_detected,
                          keytype = "COMMON",
                          columns = "ORF") %>%
      filter(!is.na(ORF)) %>%
      distinct(COMMON, ORF) %>%
      dplyr::rename(GENENAME = COMMON)
  }
)

# GO IDs for detected ORFs (org.Sc.sgd.db has GO and ONTOLOGY but not TERM)
go_ids <- AnnotationDbi::select(
  org.Sc.sgd.db,
  keys    = unique(gene_orf$ORF),
  keytype = "ORF",
  columns = c("GO", "ONTOLOGY")
) %>%
  filter(ONTOLOGY == "BP", !is.na(GO)) %>%
  dplyr::select(ORF, GO) %>%
  distinct()

# Map GO IDs → human-readable term descriptions via GO.db
go_terms <- AnnotationDbi::select(
  GO.db,
  keys    = unique(go_ids$GO),
  keytype = "GOID",
  columns = "TERM"
) %>%
  filter(!is.na(TERM))

go_bp <- go_ids %>%
  left_join(go_terms, by = c("GO" = "GOID")) %>%
  filter(!is.na(TERM)) %>%
  dplyr::select(ORF, TERM) %>%
  distinct()

# GO term → gene name sets; restrict to 3–300 genes (exclude tiny/giant terms)
go_gene_sets <- go_bp %>%
  left_join(gene_orf, by = "ORF", relationship = "many-to-many") %>%
  filter(!is.na(GENENAME)) %>%
  group_by(TERM) %>%
  summarise(genes = list(sort(unique(GENENAME))), n = n(), .groups = "drop") %>%
  filter(n >= 3, n <= 300) %>%
  arrange(TERM)

gene_universe <- list(
  genes        = genes_detected,
  go_gene_sets = go_gene_sets
)
saveRDS(gene_universe, "shiny/data/gene_universe.rds")
cat(sprintf("  gene_universe.rds: %d genes, %d GO BP terms\n",
            length(genes_detected), nrow(go_gene_sets)))

cat("\nDone. All files written to shiny/data/\n")
cat("Run the app with: shiny::runApp('shiny/')\n")
