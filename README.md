# mr.carv: An R package for performing two sample Mendelian Randomization considering both common and rare variants

## Overview

`mr.carv` is an R package designed for performing two sample Mendelian Randomization considering both common and rare variants. It leverages the functionality of the `STAARpipeline` package to provide powerful tools to get the summary statistics of whole-genome sequencing.

## Installation

You can install the released version of `mr.carv` from github with:

```
# install.packages("devtools")
devtools::install_github("yu-zhang-oYo/mr.carv")
```

## Usage

Here's a quick example to get you started:

1. load required package and the genetic association data with TC

the data are saved in the folder 'inst/extdata/' of the package

```
## load required packages
library(Matrix)
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(igraph)
library(mr.carv)

# load the indiviual variant association results
lipid_indv <- readxl::read_excel("/Users/yzh10/Library/CloudStorage/OneDrive-IndianaUniversity/research/nuMoM2b/whole_sequencing/mr.carv/inst/extdata/TC.xlsx", sheet = "indv_effect")
# load the gene centric coding region association results
lipid_gene_coding <- readxl::read_excel("/Users/yzh10/Library/CloudStorage/OneDrive-IndianaUniversity/research/nuMoM2b/whole_sequencing/mr.carv/inst/extdata/TC.xlsx", sheet = "coding")
lipid_gene_noncoding <- readxl::read_excel("/Users/yzh10/Library/CloudStorage/OneDrive-IndianaUniversity/research/nuMoM2b/whole_sequencing/mr.carv/inst/extdata/TC.xlsx", sheet = "noncoding")
lipid_window <- readxl::read_excel("/Users/yzh10/Library/CloudStorage/OneDrive-IndianaUniversity/research/nuMoM2b/whole_sequencing/mr.carv/inst/extdata/TC.xlsx", sheet = "window")
# check the data format
help(TC)

```

2. select the SNPs that are in linkage equilibrium, and save the summary statistics

Note: you can try pvalues_weight=FALSE or TRUE to check the difference of the results.

```
setwd("/Users/yzh10/Library/CloudStorage/OneDrive-IndianaUniversity/research/nuMoM2b/whole_sequencing/real_data/GA")
## Null model
obj_nullmodel <- get(load("./outputobj_nullmodel.Rdata"))

## QC_label
QC_label <- "annotation/filter"
## variant_type
variant_type_indv <- "variant"
variant_type_gene <- "SNV"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("../Annotation_name_catalog.Rdata"))
rare_maf_cutoff <- 0.05
rv_num_cutoff <- 2

results <- list()
for(chr in c(19,22)){
  ## aGDS file
  agds.path <- paste0("../gds/freeze.11a.chr", chr, ".pass_and_fail.gtonly.minDP10.gds")
  # Open the GDS file
  genofile <- seqOpen(agds.path)
  results[[paste0("chr", chr)]] <- select_LE_Estimate(chr, df_indv=lipid_indv, df_coding=lipid_gene_coding, df_noncoding=lipid_gene_noncoding, df_window=lipid_window, 
                                                      genofile=genofile, obj_nullmodel=obj_nullmodel, 
                                                      rare_maf_cutoff=rare_maf_cutoff, rv_num_cutoff=2,
                                                      QC_label=QC_label,variant_type_indv=variant_type_indv,variant_type_gene=variant_type_gene,
                                                      geno_missing_imputation=geno_missing_imputation,
                                                      Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                      le_threshold=0.1, pvalues_weight=TRUE, silent=FALSE)
  # Close the GDS file
  seqClose(genofile)
}

saveRDS(results, file = "summary_statistics.rds")

obj <- readRDS("summary_statistics.rds")
```



3. Perform the Mendelian Randomization

```
# Load necessary package
library(dplyr)
library(MendelianRandomization)

# Initialize empty data frames for the combined results
combined_X_gene <- data.frame()
combined_X_indv <- data.frame()
combined_Y_gene <- data.frame()
combined_Y_indv <- data.frame()

# Iterate over each chromosome's result
for (chromosome_result in results) {
  combined_X_gene <- bind_rows(combined_X_gene, chromosome_result$X_gene)
  combined_X_indv <- bind_rows(combined_X_indv, chromosome_result$X_indv)
  combined_Y_gene <- bind_rows(combined_Y_gene, chromosome_result$Y_gene)
  combined_Y_indv <- bind_rows(combined_Y_indv, chromosome_result$Y_indv)
}

# Create an MRInput object
MR_object <- mr_input(
  bx = c(combined_X_indv$Est, combined_X_gene$Est),
  bxse = c(combined_X_indv$Est_se, combined_X_gene$Est_se),
  by = c(combined_Y_indv$Est, combined_Y_gene$Est),
  byse = c(combined_Y_indv$Est_se, combined_Y_gene$Est_se)
)

# Perform MR analysis using IVW
MR_result <- mr_ivw(MR_object)
MR_result         
```



