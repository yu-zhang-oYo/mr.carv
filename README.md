# mr.carv: An R package for performing two sample Mendelian Randomization considering both common and rare variants

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/mr.carv)](https://cran.r-project.org/package=mr.carv)
[![Travis-CI Build Status](https://travis-ci.org/yourusername/mr.carv.svg?branch=master)](https://travis-ci.org/yourusername/mr.carv)
[![Coverage Status](https://img.shields.io/codecov/c/github/yourusername/mr.carv/master.svg)](https://codecov.io/github/yourusername/mr.carv?branch=master)

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

```
## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(igraph)
library(mr.carv)

# load the indiviual variant association results
X_indv <- readxl::read_excel("/Users/yzh10/Library/CloudStorage/OneDrive-IndianaUniversity/research/nuMoM2b/whole_sequencing/mr.carv/inst/extdata/TC_CHR19.xlsx", sheet = "indv_effect")
# load the gene centric coding region association results
X_gene <- readxl::read_excel("/Users/yzh10/Library/CloudStorage/OneDrive-IndianaUniversity/research/nuMoM2b/whole_sequencing/mr.carv/inst/extdata/TC_CHR19.xlsx", sheet = "gene_effect")

```

2. select the SNPs that are in linkage equilibrium, and save the summary statistics
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

results <- list()
for(chr in c(19,22)){
  ## aGDS file
  agds.path <- paste0("../gds/freeze.11a.chr", chr, ".pass_and_fail.gtonly.minDP10.gds")
  # Open the GDS file
  genofile <- seqOpen(agds.path)
  results[[paste0("chr", chr)]] <- select_LE_Estimate(chr, df_indv=X_indv, df_gene=X_gene, genofile=genofile, obj_nullmodel=obj_nullmodel, 
                                                      rare_maf_cutoff=rare_maf_cutoff, rv_num_cutoff=2,
                                                      QC_label=QC_label,variant_type_indv=variant_type_indv,variant_type_gene=variant_type_gene,
                                                      geno_missing_imputation=geno_missing_imputation,
                                                      Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                      le_threshold=0.1, pvalues_weight=FALSE, silent=FALSE)
  # Close the GDS file
  seqClose(genofile)
}

saveRDS(results, file = "summary_statistics.rds")

obj <- readRDS("summary_statistics.rds")
```



3. Perform the Mendelian Randomization

```

  
                 
```



