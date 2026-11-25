---
title: "01 — Data preparation and QC"
output:
  html_document:
    keep_md: true
---

**Purpose:** SILAC ratio normalization, quality control, and PCA for the
protein-level data.

**Inputs:** `data/proteinGroups.csv`

**Outputs:** `tables/master_tidy_protein.csv`, `tables/normed_ratio_perseus.csv`,
`tables/master_protein.csv`, QC box plots and PCA figures in `plots/`.

## Chunks

* [Load libraries]
* [Read in tables]
* [Proteomic data preparation]
  * [QC.1 - Intensity and normalized intensity]
    * [Data preparation]
    * [Box plot]
    * [Scatter plot]
  * [QC.2 - Ratio and normalized ratio]
    * [Data preparation]
    * [Box plot]
    * [Scatter plot]
  * [Combine ratios and intensity into one table]
* [PCA analysis for each data type]
* [T-test for DE at each timepoint]


## Load libraries

``` r
library(dplyr)
```

```
## Warning: package 'dplyr' was built under R version 4.5.2
```

``` r
library(qvalue)
```

```
## Warning: package 'qvalue' was built under R version 4.5.2
```

``` r
library(ggplot2)
```

```
## Warning: package 'ggplot2' was built under R version 4.5.2
```

``` r
library(reshape2)
library(tidyverse)
```

```
## Warning: package 'tibble' was built under R version 4.5.2
```

```
## Warning: package 'tidyr' was built under R version 4.5.2
```

```
## Warning: package 'readr' was built under R version 4.5.2
```

```
## Warning: package 'purrr' was built under R version 4.5.2
```

```
## Warning: package 'lubridate' was built under R version 4.5.2
```

``` r
library(GGally)
```

```
## Warning: package 'GGally' was built under R version 4.5.2
```

``` r
library(ggpubr)
```

```
## Warning: package 'ggpubr' was built under R version 4.5.2
```

``` r
library(ggrepel)
```

```
## Warning: package 'ggrepel' was built under R version 4.5.2
```

``` r
library(gridExtra)
```

```
## Warning: package 'gridExtra' was built under R version 4.5.2
```

``` r
library(ggfortify)
```

```
## Warning: package 'ggfortify' was built under R version 4.5.2
```

``` r
library(limma)
```

```
## Warning: package 'limma' was built under R version 4.5.2
```

``` r
library(plotly)
```

```
## Warning: package 'plotly' was built under R version 4.5.2
```

``` r
library(broom)
```

```
## Warning: package 'broom' was built under R version 4.5.2
```

``` r
library(pheatmap)
library(psych)  # geometric mean calculation
```

```
## Warning: package 'psych' was built under R version 4.5.2
```

``` r
library(RColorBrewer)
```


## Read in tables 

```
## # A tibble: 6 × 1,679
##   `Protein IDs` `Majority protein IDs` `Peptide counts (all)` `Peptide counts (razor+unique)` Peptide counts (uniqu…¹ `Protein names`
##   <chr>         <chr>                  <chr>                  <chr>                           <chr>                   <chr>          
## 1 Q9URQ5        Q9URQ5                 3                      3                               3                       High temperatu…
## 2 Q9P305        Q9P305                 11                     11                              11                      mRNA stability…
## 3 Q99394        Q99394                 4                      4                               4                       Trafficking pr…
## 4 Q99385        Q99385                 2                      2                               2                       Vacuolar calci…
## 5 Q99383        Q99383                 14                     14                              14                      Nuclear polyad…
## 6 Q99369        Q99369                 6                      6                               6                       Family of seri…
## # ℹ abbreviated name: ¹​`Peptide counts (unique)`
## # ℹ 1,673 more variables: `Gene names` <chr>, `Fasta headers` <chr>, `Number of proteins` <dbl>, Peptides <dbl>,
## #   `Razor + unique peptides` <dbl>, `Unique peptides` <dbl>, `Peptides phospho_C_T00_R1` <dbl>, `Peptides phospho_C_T00_R2` <dbl>,
## #   `Peptides phospho_C_T00_R3` <dbl>, `Peptides phospho_C_T06_R1` <dbl>, `Peptides phospho_C_T06_R2` <dbl>,
## #   `Peptides phospho_C_T06_R3` <dbl>, `Peptides phospho_C_T16_R1` <dbl>, `Peptides phospho_C_T16_R2` <dbl>,
## #   `Peptides phospho_C_T16_R3` <dbl>, `Peptides phospho_C_T30_R1` <dbl>, `Peptides phospho_C_T30_R2` <dbl>,
## #   `Peptides phospho_C_T30_R3` <dbl>, `Peptides phospho_N_T06_R1` <dbl>, `Peptides phospho_N_T06_R2` <dbl>, …
```

```
## [1] 3452 1679
```


## Proteomic data preparation
* This code extract the information from proteinGroup file generated in MaxQuant and put them in a tidy format
  * intensity
    * Run Anders & Huber, 2010 method to get normalized intensity for each protein
  * ratio
  * normalized_ratio

### QC.1 - Intensity and normalized intensity
#### Data preparation


#### Box plot
![](01_data_preparation_qc_files/figure-html/unnamed-chunk-4-1.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-4-2.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-4-3.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-4-4.png)<!-- -->

#### Scatter plot
![](01_data_preparation_qc_files/figure-html/unnamed-chunk-5-1.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-5-2.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-5-3.png)<!-- -->

### QC.2 - Ratio and normalized ratio
#### Data preparation


#### Box plot
![](01_data_preparation_qc_files/figure-html/unnamed-chunk-7-1.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-7-2.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-7-3.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-7-4.png)<!-- -->

#### Scatter plot
![](01_data_preparation_qc_files/figure-html/unnamed-chunk-8-1.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-8-2.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-8-3.png)<!-- -->

### Combine ratios and intensity into one table 

```
## [1] 551664     10
```

```
##   Protein IDs Gene names                                 Protein names
## 1      Q9P305       IGO2                   mRNA stability protein IGO2
## 2      Q99385       VCX1              Vacuolar calcium ion transporter
## 3      Q99383       HRP1  Nuclear polyadenylated RNA-binding protein 4
## 4      Q99297       ODC2    Mitochondrial 2-oxodicarboxylate carrier 2
## 5      Q99258       RIB3 3,4-dihydroxy-2-butanone 4-phosphate synthase
## 6      Q99257      MEX67                      mRNA export factor MEX67
##                                                                                                                                             Fasta headers
## 1                   sp|Q9P305|IGO2_YEAST mRNA stability protein IGO2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=IGO2 PE=1 SV=3
## 2              sp|Q99385|VCX1_YEAST Vacuolar calcium ion transporter OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=VCX1 PE=1 SV=1
## 3  sp|Q99383|HRP1_YEAST Nuclear polyadenylated RNA-binding protein 4 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=HRP1 PE=1 SV=1
## 4    sp|Q99297|ODC2_YEAST Mitochondrial 2-oxodicarboxylate carrier 2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=ODC2 PE=1 SV=1
## 5 sp|Q99258|RIB3_YEAST 3,4-dihydroxy-2-butanone 4-phosphate synthase OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=RIB3 PE=1 SV=1
## 6                    sp|Q99257|MEX67_YEAST mRNA export factor MEX67 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=MEX67 PE=1 SV=1
##   nutrient samplingTime replicate genotype      data     value
## 1        C          T00        R1 spike-in intensity  35919000
## 2        C          T00        R1 spike-in intensity  27087000
## 3        C          T00        R1 spike-in intensity 306210000
## 4        C          T00        R1 spike-in intensity  82491000
## 5        C          T00        R1 spike-in intensity 157850000
## 6        C          T00        R1 spike-in intensity  85041000
```

## PCA analysis for each data type
![](01_data_preparation_qc_files/figure-html/unnamed-chunk-10-1.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-10-2.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-10-3.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-10-4.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-10-5.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-10-6.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-10-7.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-10-8.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-10-9.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-10-10.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-10-11.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-10-12.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-10-13.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-10-14.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-10-15.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-10-16.png)<!-- -->

#### clearly, normalized ratio maximized the difference on proteome proflies between genotypes

## Filter out the poorly correlated replicates
#### 1. WT -C T06 rep1, 
#### 2. WT -C T16 rep1, 
#### 3. WT -C T30 rep2, 
#### 4. rim15 -C T16 rep2, 
#### 5. rim15 -C T30 rep2, 
#### 6. WT -P T16 rep1, 
#### 7. rim15KO -P T16 rep1

``` r
master <- read_csv("./tables/master_protein.csv")
```

```
## Rows: 551664 Columns: 10
## ── Column specification ─────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Delimiter: ","
## chr (9): Protein IDs, Gene names, Protein names, Fasta headers, nutrient, samplingTime, replicate, genotype, data
## dbl (1): value
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

``` r
filtered <- master %>% 
    #remove spike-in in the data frame
  filter(genotype != "spike-in") %>%
    #remove samples from nitrogen starved condition
  filter(nutrient != "N") %>%
    #filter(replicate != "R2") %>%
  filter(data == "normed_ratio") %>% 
  melt() %>%
  dplyr::select(-variable) %>%
  #remove na values
  dplyr::filter(!is.na(value)) %>% 
  #remove infinite values
  dplyr::filter(is.finite(value)) %>% 
  #remove 0 values, wich will be a problem for later log based transformation
  dplyr::filter(value > 0) %>%
  group_by(`Protein IDs`,`Gene names`, nutrient, samplingTime, replicate, data) %>%
  mutate(genotype_number = n()) %>%
  filter(genotype_number == 2) %>%
  group_by(`Protein IDs`,`Gene names`, nutrient, samplingTime, genotype, data) %>%
  mutate(rep_numner = n()) %>%
  filter(rep_numner == 3) %>%
  group_by(`Protein IDs`,`Gene names`, nutrient, replicate, genotype, data) %>%
  mutate(time_number = n()) %>%
  filter(time_number == 4) %>%
  filter(!(
    #1. WT -C T06 rep1, 
    (genotype == 'WT' & nutrient == "C" & samplingTime == 'T06' & replicate == 'R1') |
      #WT -C T16 rep1
      (genotype == 'WT' & nutrient == "C" & samplingTime == 'T16' & replicate == 'R1') |
      #WT -C T30 rep2, 
      (genotype == 'WT' & nutrient == "C" & samplingTime == 'T30' & replicate == 'R2') |
      #rim15 -C T16 rep2, 
      (genotype == 'rim15KO' & nutrient == "C" & samplingTime == 'T16' & replicate == 'R2') |
      #rim15 -C T30 rep2, 
      (genotype == 'rim15KO' & nutrient == "C" & samplingTime == '30' & replicate == 'R2') |
      #6. WT -P T16 rep1, 
      (genotype == 'WT' & nutrient == "P" & samplingTime == 'T16' & replicate == 'R1') |
      #7. rim15KO -P T16 rep1
      (genotype == 'rim15KO' & nutrient == "P" & samplingTime == 'T16' & replicate == 'R1') ))
```

```
## Using Protein IDs, Gene names, Protein names, Fasta headers, nutrient, samplingTime, replicate, genotype, data as id variables
```

``` r
#save
write_csv (filtered, "./tables/filtered_norm_ratios_WP.csv")
```

## T-test for DE at each timepoint 
![](01_data_preparation_qc_files/figure-html/unnamed-chunk-12-1.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-12-2.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-12-3.png)<!-- -->![](01_data_preparation_qc_files/figure-html/unnamed-chunk-12-4.png)<!-- -->
