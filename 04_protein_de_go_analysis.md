---
title: "04 — Differential protein abundance and GO/KEGG enrichment"
output:
  html_document:
    keep_md: true
---

**Purpose:** ANCOVA-based differential expression by genotype and nutrient;
GO and KEGG enrichment of significant proteins; visualization (volcano plots,
heatmaps, enrichment bar charts).

**Inputs:** `tables/master_tidy_protein.csv`, `tables/normed_ratio_perseus.csv`

**Outputs:** `tables/DE_proteins_byGenotype.csv`, `tables/DE_proteins_byNutrient.csv`,
`tables/anova_genotype.csv`, `tables/anova_nutrient.csv`,
`tables/anova_genotype_nutrient.csv`, figures in `plots/`.

## Chrunks

  * [Load libraries]
  * [Read in pre-processed protein table]
* [Two-way ANCOVA]
  * [Genotype wise]
    * [Modeling]
    * [Summary of model fitting]
    * [DE visualization]
    * [Kegg analysis visualization]
  * [Nutrient wise]
    * [Modeling]
    * [Summary of model fitting]
    * [DE visualization]
    * [Kegg analysis visualization]
* [Three way ANCOVA]
  * [Modeling]
  * [Summary of model fitting]
  * [DE visualization]
  * [Kegg analysis visualization]
* [Single protein check]

## Load libraries

 
## Read in pre-processed protein table


# Two-way ANCOVA
## Genotype wise
### Modeling


### Summary of model fitting

```
## Rows: 34304 Columns: 10
## ── Column specification ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Delimiter: ","
## chr (5): nutrient, Gene names, Protein IDs, data, term
## dbl (5): estimate, std.error, statistic, p.value, qvalue
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
## `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```
## `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-4-2.png)<!-- -->

```
## `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-4-3.png)<!-- -->

```
## [1] 384  11
```

### DE visualization

``` r
summary <- sig_genotype %>% 
  filter(data %in% c("normed_intensity","normed_ratio")) %>%
  group_by(nutrient, data, type) %>%
  summarize(DE = n())
```

```
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by nutrient, data, and type.
## ℹ Output is grouped by nutrient and data.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(nutrient, data, type))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
```

``` r
print(summary)
```

```
## # A tibble: 7 × 4
## # Groups:   nutrient, data [4]
##   nutrient data             type                    DE
##   <chr>    <chr>            <chr>                <int>
## 1 C        normed_intensity upregulated in wt        6
## 2 C        normed_ratio     upregulated in rim15   102
## 3 C        normed_ratio     upregulated in wt       23
## 4 P        normed_intensity upregulated in rim15     6
## 5 P        normed_intensity upregulated in wt       10
## 6 P        normed_ratio     upregulated in rim15    34
## 7 P        normed_ratio     upregulated in wt       21
```

``` r
listInput <- list(
  `upregulated in wt (-C)` = sig_genotype %>% filter(nutrient == "C" & data == "normed_ratio" & estimate > 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(), 
  `upregulated in wt (-P)` = sig_genotype %>% filter(nutrient == "P" & data == "normed_ratio" & estimate > 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(),
  `upregulated in rim15 (-C)` = sig_genotype %>% filter(nutrient == "C" & data == "normed_ratio" & estimate < 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(),
  `upregulated in rim15 (-P)` = sig_genotype %>% filter(nutrient == "P" & data == "normed_ratio" & estimate < 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector())

library("SuperExactTest")
set.seed(1234)
n=1074

sig_geno = supertest(listInput, n=n)

summary <- summary(sig_geno)$Table 

ol_genes <- summary %>% filter(P.value < 0.05) %>% dplyr::select(Elements)
  
plot(sig_geno, Layout="landscape", sort.by="size", keep=FALSE,
	bar.split=c(70,180), elements.cex=0.7, #show.elements=TRUE, 
	elements.list=subset(summary(sig_geno)$Table, Observed.Overlap <= 200),
	show.expected.overlap=TRUE, expected.overlap.style="hatchedBox",	color.expected.overlap='#a6beff')
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

``` r
upset(fromList(listInput), order.by = "freq")
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-5-2.png)<!-- -->

``` r
###normalized intensity
listInput <- list(
  `upregulated in wt (-C)` = sig_genotype %>% filter(nutrient == "C" & data == "normed_intensity" & estimate > 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(), 
  `upregulated in wt (-P)` = sig_genotype %>% filter(nutrient == "P" & data == "normed_intensity" & estimate > 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(),
  `upregulated in rim15 (-C)` = sig_genotype %>% filter(nutrient == "C" & data == "normed_intensity" & estimate < 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(),
  `upregulated in rim15 (-P)` = sig_genotype %>% filter(nutrient == "P" & data == "normed_intensity" & estimate < 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector())

sig_geno = supertest(listInput, n=n)

summary <- summary(sig_geno)$Table 

ol_genes <- summary %>% filter(P.value < 0.05) %>% dplyr::select(Elements)
  
plot(sig_geno, Layout="landscape", sort.by="size", keep=FALSE,
	bar.split=c(70,180), elements.cex=0.7, #show.elements=TRUE, 
	elements.list=subset(summary(sig_geno)$Table, Observed.Overlap <= 200),
	show.expected.overlap=TRUE, expected.overlap.style="hatchedBox",	color.expected.overlap='#a6beff')
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-5-3.png)<!-- -->

``` r
upset(fromList(listInput), order.by = "freq")
```

```
## `geom_line()`: Each group consists of only one observation.
## ℹ Do you need to adjust the group aesthetic?
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-5-4.png)<!-- -->

### Kegg analysis visualization

``` r
sig_genotype <- read_csv("./tables/anova_genotype_sig.csv")
```

```
## Rows: 384 Columns: 11
## ── Column specification ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Delimiter: ","
## chr (6): nutrient, Gene names, Protein IDs, data, term, type
## dbl (5): estimate, std.error, statistic, p.value, qvalue
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

``` r
#convert protein IDs column into ENSEMBL, which is readable for kegg annotation
tmp<-bitr(sig_genotype$`Protein IDs`, "UNIPROT", "ENSEMBL", "org.Sc.sgd.db", drop = TRUE)
```

```
## 'select()' returned 1:1 mapping between keys and columns
```

```
## Warning in bitr(sig_genotype$`Protein IDs`, "UNIPROT", "ENSEMBL", "org.Sc.sgd.db", : 3.38% of input gene IDs are fail to map...
```

``` r
sig <- left_join (tmp, sig_genotype %>% dplyr::rename(UNIPROT = `Protein IDs`), by = "UNIPROT")

sig_C_upWT <- sig %>% 
  filter(data == "normed_ratio" & estimate > 0) %>% 
  dplyr::select(nutrient, ENSEMBL, UNIPROT)

upWT <- compareCluster(ENSEMBL ~ nutrient, 
                           data=sig_C_upWT, 
                           organism = "sce",
                           fun="enrichKEGG",
                           pAdjustMethod = "BH",
                           qvalueCutoff=0.05)
```

```
## Warning in compareCluster(ENSEMBL ~ nutrient, data = sig_C_upWT, organism = "sce", : No enrichment found in any of gene cluster, please
## check your input...
```

``` r
if (!is.null(upWT)) {
  print(dotplot(upWT) +
    ggplot2::theme(axis.text.x=element_text(angle=90),legend.position="bottom"))
} else {
  message("upWT: no KEGG enrichment results (compareCluster returned NULL)")
}
```

```
## upWT: no KEGG enrichment results (compareCluster returned NULL)
```

``` r
#ggsave("./plots/go_up_C_kegg.pdf", width=10,height=8)

sig_upRim15 <- sig %>% 
  filter(data == "normed_ratio" & estimate < 0) %>% 
  dplyr::select(nutrient, ENSEMBL, UNIPROT)

upRim15 <- compareCluster(ENSEMBL ~ nutrient, 
                           data=sig_upRim15,
                           organism = "sce",
                           #OrgDb=org.Sc.sgd.db::org.Sc.sgd.db,
                           #ont="BP", 
                           fun="enrichKEGG", 
                           #keyType = "UNIPROT", 
                           pAdjustMethod = "BH",
                           qvalueCutoff=0.05)

if (!is.null(upRim15)) {
  print(dotplot(upRim15) +
    ggplot2::theme(axis.text.x=element_text(angle=90),legend.position="bottom"))
} else {
  message("upRim15: no KEGG enrichment results (compareCluster returned NULL)")
}
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

``` r
#ggsave("./plots/go_up_P.pdf", width=10,height=8)
```

## Nutrient wise
### Modeling


### Summary of model fitting

```
## Rows: 34304 Columns: 10
## ── Column specification ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Delimiter: ","
## chr (5): genotype, Gene names, Protein IDs, data, term
## dbl (5): estimate, std.error, statistic, p.value, qvalue
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
## `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```
## `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-8-2.png)<!-- -->

```
## `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-8-3.png)<!-- -->

```
## [1] 2300   11
```

```
## # A tibble: 6 × 11
##   genotype `Gene names` `Protein IDs` data             term                   estimate std.error statistic p.value  qvalue type          
##   <chr>    <chr>        <chr>         <chr>            <chr>                     <dbl>     <dbl>     <dbl>   <dbl>   <dbl> <chr>         
## 1 rim15KO  AAP1         P37898        normed_intensity samplingTime:nutrientP  -0.0480   0.0131      -3.67 0.00151 0.00337 upregulated i…
## 2 rim15KO  AAP1         P37898        normed_ratio     samplingTime:nutrientP  -0.0231   0.00778     -2.97 0.00753 0.0139  upregulated i…
## 3 rim15KO  AAP1         P37898        ratio            samplingTime:nutrientP  -0.0229   0.00938     -2.44 0.0240  0.0381  upregulated i…
## 4 rim15KO  AAT2         P23542        normed_intensity samplingTime:nutrientP  -0.0388   0.0152      -2.55 0.0189  0.0310  upregulated i…
## 5 rim15KO  AAT2         P23542        normed_ratio     samplingTime:nutrientP  -0.0468   0.0124      -3.77 0.00122 0.00277 upregulated i…
## 6 rim15KO  AAT2         P23542        ratio            samplingTime:nutrientP  -0.0362   0.0147      -2.46 0.0232  0.0369  upregulated i…
```

### DE visualization
#### normalized ratio

``` r
summary <- sig_nutrient %>% 
  filter(data %in% c("normed_intensity","normed_ratio")) %>%
  group_by(genotype, data, type) %>%
  summarize(DE = n())
```

```
## `summarise()` has regrouped the output.
## ℹ Summaries were computed grouped by genotype, data, and type.
## ℹ Output is grouped by genotype and data.
## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
## ℹ Use `summarise(.by = c(genotype, data, type))` for per-operation grouping (`?dplyr::dplyr_by`) instead.
```

``` r
print(summary)
```

```
## # A tibble: 8 × 4
## # Groups:   genotype, data [4]
##   genotype data             type                DE
##   <chr>    <chr>            <chr>            <int>
## 1 WT       normed_intensity upregulated in C    78
## 2 WT       normed_intensity upregulated in P   242
## 3 WT       normed_ratio     upregulated in C     8
## 4 WT       normed_ratio     upregulated in P   415
## 5 rim15KO  normed_intensity upregulated in C   108
## 6 rim15KO  normed_intensity upregulated in P   268
## 7 rim15KO  normed_ratio     upregulated in C    13
## 8 rim15KO  normed_ratio     upregulated in P   512
```

``` r
listInput <- list(
  `upregulated in C (wt)` = sig_nutrient %>% filter(genotype == "WT" & data == "normed_ratio" & estimate > 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(), 
  `upregulated in P (wt)` = sig_nutrient %>% filter(genotype == "WT" & data == "normed_ratio" & estimate < 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(),
  `upregulated in C (rim15)` = sig_nutrient %>% filter(genotype == "rim15KO" & data == "normed_ratio" & estimate > 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(),
  `upregulated in P (rim15)` = sig_nutrient %>% filter(genotype == "rim15KO" & data == "normed_ratio" & estimate < 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector())

library("SuperExactTest")
set.seed(1234)
n=1074

sig = supertest(listInput, n=n)

summary <- summary(sig)$Table 

ol_genes <- summary %>% filter(P.value < 0.05)
#%>% select(Elements)
  
plot(sig, Layout="landscape", sort.by="size", keep=FALSE,	#bar.split=c(70,180), 
     elements.cex=0.7, #show.elements=TRUE, 
	elements.list=subset(summary(sig)$Table, Observed.Overlap <= 200),
	show.expected.overlap=TRUE, expected.overlap.style="hatchedBox",	color.expected.overlap='#a6beff')
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

``` r
upset(fromList(listInput), order.by = "freq")
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-9-2.png)<!-- -->

#### normalized intensity

``` r
listInput <- list(
  `upregulated in C (wt)` = sig_nutrient %>% filter(genotype == "WT" & data == "normed_intensity" & estimate > 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(), 
  `upregulated in P (wt)` = sig_nutrient %>% filter(genotype == "WT" & data == "normed_intensity" & estimate < 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(),
  `upregulated in C (rim15)` = sig_nutrient %>% filter(genotype == "rim15KO" & data == "normed_intensity" & estimate > 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector(),
  `upregulated in P (rim15)` = sig_nutrient %>% filter(genotype == "rim15KO" & data == "normed_intensity" & estimate < 0) %>% dplyr::select(`Protein IDs`) %>% unlist() %>% as.vector())

sig = supertest(listInput, n=n)

summary <- summary(sig)$Table 

ol_genes <- summary %>% filter(P.value < 0.05)
#%>% select(Elements)
  
plot(sig, Layout="landscape", sort.by="size", keep=FALSE,	#bar.split=c(70,180), 
     elements.cex=0.7, #show.elements=TRUE, 
	elements.list=subset(summary(sig)$Table, Observed.Overlap <= 200),
	show.expected.overlap=TRUE, expected.overlap.style="hatchedBox",	color.expected.overlap='#a6beff')
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

``` r
upset(fromList(listInput), order.by = "freq")
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-10-2.png)<!-- -->

### Kegg analysis & visualization

``` r
sig_nutrient <- read_csv("./tables/ancova_nutrient_sig.csv")
```

```
## Rows: 2300 Columns: 11
## ── Column specification ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Delimiter: ","
## chr (6): genotype, Gene names, Protein IDs, data, term, type
## dbl (5): estimate, std.error, statistic, p.value, qvalue
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

``` r
#convert protein IDs column into ENSEMBL, which is readable for kegg annotation
tmp<-bitr(sig_nutrient$`Protein IDs`, "UNIPROT", "ENSEMBL", "org.Sc.sgd.db", drop = TRUE)
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```
## Warning in bitr(sig_nutrient$`Protein IDs`, "UNIPROT", "ENSEMBL", "org.Sc.sgd.db", : 5.36% of input gene IDs are fail to map...
```

``` r
sig <- left_join (tmp, sig_nutrient %>% dplyr::rename(UNIPROT = `Protein IDs`), by = "UNIPROT")
```

```
## Warning in left_join(tmp, sig_nutrient %>% dplyr::rename(UNIPROT = `Protein IDs`), : Detected an unexpected many-to-many relationship between `x` and `y`.
## ℹ Row 1 of `x` matches multiple rows in `y`.
## ℹ Row 238 of `y` matches multiple rows in `x`.
## ℹ If a many-to-many relationship is expected, set `relationship = "many-to-many"` to silence this warning.
```

``` r
sig_nutrient_genelist_upC <- sig %>% 
  filter(data == "normed_ratio" & estimate > 0) %>% 
  dplyr::select(genotype, ENSEMBL, UNIPROT)

wt_rim15_upC <- compareCluster(ENSEMBL ~ genotype, 
                           data=sig_nutrient_genelist_upC, 
                           organism = "sce",
                           fun="enrichKEGG",
                           pAdjustMethod = "BH",
                           qvalueCutoff=0.05)

if (!is.null(wt_rim15_upC)) {
  print(dotplot(wt_rim15_upC) +
    ggplot2::theme(axis.text.x=element_text(angle=90),legend.position="bottom"))
} else {
  message("wt_rim15_upC: no KEGG enrichment results (compareCluster returned NULL)")
}
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

``` r
#ggsave("./plots/go_up_C_kegg.pdf", width=10,height=8)

sig_nutrient_genelist_upP <- sig %>% 
  filter(data == "normed_ratio" & estimate < 0) %>% 
  dplyr::select(genotype, ENSEMBL, UNIPROT)

wt_rim15_upP <- compareCluster(ENSEMBL ~ genotype, 
                           data=sig_nutrient_genelist_upP,
                           organism = "sce",
                           #OrgDb=org.Sc.sgd.db::org.Sc.sgd.db,
                           #ont="BP", 
                           fun="enrichKEGG", 
                           #keyType = "UNIPROT", 
                           pAdjustMethod = "BH",
                           qvalueCutoff=0.05)

if (!is.null(wt_rim15_upP)) {
  print(dotplot(wt_rim15_upP) +
    ggplot2::theme(axis.text.x=element_text(angle=90),legend.position="bottom"))
} else {
  message("wt_rim15_upP: no KEGG enrichment results (compareCluster returned NULL)")
}
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-11-2.png)<!-- -->

``` r
#ggsave("./plots/go_up_P.pdf", width=10,height=8)
```


## Three way ANCOVA
### Modeling


### Summary of model fitting

```
## Rows: 30016 Columns: 10
## ── Column specification ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Delimiter: ","
## chr (4): Protein IDs, Gene names, data, term
## dbl (6): df, sumsq, meansq, statistic, p.value, qvalue
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
## `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

```
## `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-13-2.png)<!-- -->

```
## `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-13-3.png)<!-- -->

```
## [1] 86 10
```

```
## # A tibble: 6 × 10
##   `Protein IDs` `Gene names` data         term                              df  sumsq meansq statistic  p.value  qvalue
##   <chr>         <chr>        <chr>        <chr>                          <dbl>  <dbl>  <dbl>     <dbl>    <dbl>   <dbl>
## 1 O14467        MBF1         normed_ratio samplingTime:genotype:nutrient     1 0.165  0.165       7.59 0.00879  0.0248 
## 2 P00445        SOD1         normed_ratio samplingTime:genotype:nutrient     1 0.357  0.357      13.9  0.000593 0.00224
## 3 P00447        SOD2         normed_ratio samplingTime:genotype:nutrient     1 0.940  0.940       7.49 0.00922  0.0259 
## 4 P00447        SOD2         ratio        samplingTime:genotype:nutrient     1 1.27   1.27        8.32 0.00628  0.0185 
## 5 P01120        RAS2         normed_ratio samplingTime:genotype:nutrient     1 0.271  0.271       7.09 0.0111   0.0303 
## 6 P04802        DPS1         normed_ratio samplingTime:genotype:nutrient     1 0.0851 0.0851      8.25 0.00648  0.0190
```

### DE expression visualization

``` r
sig_3way <- read_csv("./tables/ancova_nutrient_sig.csv")
```

```
## Rows: 2300 Columns: 11
## ── Column specification ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Delimiter: ","
## chr (6): genotype, Gene names, Protein IDs, data, term, type
## dbl (5): estimate, std.error, statistic, p.value, qvalue
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

``` r
sig_3way_proteins <- sig_3way %>% filter(data == "normed_ratio")
gene_list <- sig_3way_proteins$`Gene names`

nr_df <- master %>% 
    #remove spike-in in the data frame
    filter(genotype != "spike-in") %>%
    #remove samples from nitrogen starved condition
    filter(nutrient != "N") %>%
    #filter(replicate != "R2") %>%
    filter(data == "normed_ratio") %>%
    dplyr::select(`Protein IDs`, `Gene names`, genotype, nutrient, replicate, samplingTime, value) %>%
    filter(!is.na(value)) %>%
    filter(is.finite(value)) %>%
  filter(`Gene names` %in% sig_3way_proteins$`Gene names`) %>%
    #filter(value > 0 ) %>%
    unite(cond, nutrient, genotype, samplingTime, replicate) %>%
    spread(key = cond, value = value) 

annotation_col <- names(nr_df)[-c(1:2)] %>% as.tibble() %>%
    separate(value, c("nutrient","genotype","samplingTime","rep"), sep = "_", remove = FALSE) %>%
    column_to_rownames(var = "value")
  
  breaksList = seq(-0.07, 0.07, by = 0.0001)
#make heatmaps
  pheatmap(nr_df[,-c(1:2)] %>% drop_na(), 
           colorRampPalette(rev(brewer.pal(n =10, name = "RdYlBu")))(length(breaksList)),
           cluster_cols = T, scale = "row", 
           annotation_col = annotation_col)
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

### Kegg analysis visualization

``` r
#GO term analysis
ego_3way <- clusterProfiler::enrichGO(gene = sig_3way_proteins$`Protein IDs`, 
             ont ="ALL", 
             universe = unique(master$`Protein IDs`),
             keyType = "UNIPROT", 
             OrgDb = "org.Sc.sgd.db", 
             pAdjustMethod = "BH", 
             qvalueCutoff = 0.05)
 
if (!is.null(ego_3way)) {
  print(dotplot(ego_3way, ~ONTOLOGY))
  dir.create("./tables/GOterm", showWarnings = FALSE, recursive = TRUE)
  write_csv(as.data.frame(ego_3way),"./tables/GOterm/ego_3way_ANCOVA.csv")
} else {
  message("ego_3way: no GO enrichment results (enrichGO returned NULL)")
}
```

![](04_protein_de_go_analysis_files/figure-html/unnamed-chunk-15-1.png)<!-- -->
# Single protein check 
## For DE proteins from three-way ancova


