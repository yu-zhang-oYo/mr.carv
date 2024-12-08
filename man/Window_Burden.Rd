% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Window_Burden.R
\name{Window_Burden}
\alias{Window_Burden}
\title{The \code{Window_Burden} function performs burden analysis for each sliding/dynamic window region in the given dataframe.}
\usage{
Window_Burden(
  chr,
  df_window,
  genofile,
  obj_nullmodel,
  rare_maf_cutoff = 0.01,
  rv_num_cutoff = 2,
  QC_label = "annotation/filter",
  variant_type = c("SNV", "Indel", "variant"),
  geno_missing_imputation = c("mean", "minor"),
  Annotation_dir = "annotation/info/FunctionalAnnotation",
  Annotation_name_catalog,
  silent = FALSE
)
}
\arguments{
\item{chr}{Chromosome number.}

\item{df_window}{Dataframe containing set information and sets' association result with the exposure. Columns include CHR, Start.Loc, End.Loc, annotation, Est, Est_se, pvalue.}

\item{genofile}{An object of an opened annotated GDS (aGDS) file.}

\item{obj_nullmodel}{An object from fitting the null model, which is either the output from the \code{\link{fit_nullmodel}} function,
or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{\link{genesis2staar_nullmodel}} function.}

\item{rare_maf_cutoff}{Minor allele frequency cutoff for rare variants. Default is 0.01.}

\item{rv_num_cutoff}{Cutoff for the minimum number of variants in a given variant-set. Default is 2.}

\item{QC_label}{Channel name of the QC label in the GDS/aGDS file. Default is "annotation/filter".}

\item{variant_type}{Type of variants to include in the analysis. Choices include "SNV", "Indel", or "variant". Default is c("SNV", "Indel", "variant").}

\item{geno_missing_imputation}{Method for imputing missing genotypes. Choices are "mean" or "minor". Default is c("mean", "minor").}

\item{Annotation_dir}{Directory for functional annotation information in the aGDS file. Default is "annotation/info/FunctionalAnnotation".}

\item{Annotation_name_catalog}{A data frame containing the name and the corresponding channel name in the aGDS file.}

\item{silent}{Logical; if TRUE, suppresses messages. Default is FALSE.}
}
\value{
A list of results from sliding/dynamic window burden analysis for each set.
}
\description{
The \code{Window_Burden} function performs burden analysis for each sliding/dynamic window region in the given dataframe.
}
