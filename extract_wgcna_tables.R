# extract_wgcna_tables.R
# Run from the repo root:  Rscript extract_wgcna_tables.R
#
# Reads the four WGCNA RData files in output/wgcna/ and writes two tables:
#
#   tables/wgcna_module_membership.csv
#     One row per protein × network combination.  Columns: protein_id, gene_name,
#     nutrient, genotype, module.  Grey (WGCNA's unassigned bin) is labelled
#     "unassigned".  A wide version (one column per network) is also written:
#     tables/wgcna_module_membership_wide.csv
#
#   tables/wgcna_module_eigengenes.csv
#     Module eigengene values per sample per network, long format.  Columns:
#     network, nutrient, genotype, timepoint_h, replicate, module, eigengene.
#     Grey/unassigned eigengenes are excluded.

suppressPackageStartupMessages(library(tidyverse))

# ── paths ────────────────────────────────────────────────────────────────────
rdata_files <- c(
  C_WT    = "output/wgcna/C_WT_ratio_networkConstruction-auto.RData",
  C_rim15 = "output/wgcna/C_rim15KO_ratio_networkConstruction-auto.RData",
  P_WT    = "output/wgcna/P_WT_ratio_networkConstruction-auto.RData",
  P_rim15 = "output/wgcna/P_rim15KO_ratio_networkConstruction-auto.RData"
)

stopifnot(all(file.exists(rdata_files)))

# Gene-name lookup from master_protein.csv (Protein IDs → first gene name)
gene_map <- read_csv("tables/master_protein.csv", show_col_types = FALSE) %>%
  distinct(`Protein IDs`, `Gene names`) %>%
  rename(protein_id = `Protein IDs`, gene_name = `Gene names`) %>%
  mutate(gene_name = str_extract(gene_name, "^[^;]+"))   # first name in group

# ── helper: numeric module label → color name ─────────────────────────────
# WGCNA::labels2colors maps 0→grey, 1→turquoise, 2→blue, 3→brown, …
# Replicate the mapping without requiring the WGCNA package at extraction time.
label_to_color <- function(labels) {
  palette <- c("grey", "turquoise", "blue", "brown", "yellow", "green",
               "red", "black", "pink", "magenta", "purple", "greenyellow",
               "tan", "salmon", "cyan", "midnightblue", "lightcyan",
               "grey60", "lightgreen", "lightyellow", "royalblue",
               "darkred", "darkgreen", "darkturquoise", "darkgrey",
               "orange", "darkorange", "white", "skyblue", "saddlebrown",
               "steelblue", "paleturquoise", "violet", "darkolivegreen",
               "darkmagenta")
  palette[labels + 1]   # labels are 0-indexed; +1 for R 1-indexing
}

# ── 1. Module membership (long) ───────────────────────────────────────────
cat("Extracting module membership...\n")

membership_long <- map_dfr(names(rdata_files), function(nm) {
  e <- new.env()
  load(rdata_files[nm], envir = e)

  # moduleLabels is a named integer vector: name = Protein ID, value = module number
  tibble(
    network    = nm,
    nutrient   = str_extract(nm, "^[CP]"),
    genotype   = str_remove(nm, "^[CP]_"),
    protein_id = names(e$moduleLabels),
    module_num = as.integer(e$moduleLabels),
    module     = label_to_color(as.integer(e$moduleLabels))
  ) %>%
    mutate(module = if_else(module == "grey", "unassigned", module))
}) %>%
  left_join(gene_map, by = "protein_id") %>%
  select(network, nutrient, genotype, protein_id, gene_name, module)

write_csv(membership_long, "tables/wgcna_module_membership.csv")
cat(sprintf("  wgcna_module_membership.csv: %d rows\n", nrow(membership_long)))

# Module sizes (excluding unassigned)
cat("\nModule sizes (assigned proteins only):\n")
membership_long %>%
  filter(module != "unassigned") %>%
  count(network, module) %>%
  print(n = 20)

# ── 2. Wide membership table: one column per network ─────────────────────
cat("\nBuilding wide membership table...\n")

membership_wide <- membership_long %>%
  select(protein_id, gene_name, network, module) %>%
  pivot_wider(names_from = network, values_from = module,
              names_glue = "{network}_module") %>%
  arrange(protein_id)

write_csv(membership_wide, "tables/wgcna_module_membership_wide.csv")
cat(sprintf("  wgcna_module_membership_wide.csv: %d proteins × %d columns\n",
            nrow(membership_wide), ncol(membership_wide)))

# ── 3. Module eigengenes (long) ───────────────────────────────────────────
cat("\nExtracting module eigengenes...\n")

eigengenes_long <- map_dfr(names(rdata_files), function(nm) {
  e <- new.env()
  load(rdata_files[nm], envir = e)

  # MEs: rows = samples ("C_WT_T00_R1"), cols = "ME0", "ME1", ...
  # Map numeric ME columns → color names, drop grey (ME0)
  me_df <- as.data.frame(e$MEs) %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "me_label", values_to = "eigengene") %>%
    mutate(
      module_num = as.integer(str_remove(me_label, "^ME")),
      module     = label_to_color(module_num)
    ) %>%
    filter(module != "grey") %>%    # drop unassigned eigengene
    separate(sample, into = c("nutrient_s", "genotype_s", "timepoint", "replicate"),
             sep = "_") %>%
    mutate(
      timepoint_h = as.integer(str_remove(timepoint, "^T")),
      network     = nm,
      nutrient    = str_extract(nm, "^[CP]"),
      genotype    = str_remove(nm, "^[CP]_")
    ) %>%
    select(network, nutrient, genotype, timepoint_h, replicate, module, eigengene)
}) %>%
  arrange(network, module, timepoint_h, replicate)

write_csv(eigengenes_long, "tables/wgcna_module_eigengenes.csv")
cat(sprintf("  wgcna_module_eigengenes.csv: %d rows\n", nrow(eigengenes_long)))

cat("\nDone. Files written to tables/:\n")
cat("  wgcna_module_membership.csv        (long: one row per protein × network)\n")
cat("  wgcna_module_membership_wide.csv   (wide: one row per protein, columns = networks)\n")
cat("  wgcna_module_eigengenes.csv        (long: eigengene per sample × module × network)\n")
