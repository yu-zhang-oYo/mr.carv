% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TC}
\alias{TC}
\title{TC Data}
\format{
A list with four data frames:
\describe{
  \item{indv_effect}{A data frame with columns:
    \itemize{
      \item \code{CHR}: Numeric, chromosome
      \item \code{POS}: Numeric, the position of each individual variant
      \item \code{REF}: Character, reference allele of each individual variant
      \item \code{ALT}: Character, alternative allele of each individual variant
      \item \code{annotation.id}: Character, the annotation id of each individual variant
      \item \code{Est}: Numeric, the effect size estimates of each individual variant on the exposure total cholesterol
      \item \code{Est_se}: Numeric, the standard error of the effect size estimates for each individual variant on the exposure total cholesterol
      \item \code{pvalue}: Numeric, the P value of the effect size for each individual variant
      \item \code{MAF}: Numeric, optional, for the harmonization
    }
  }
  \item{coding}{A data frame with columns:
    \itemize{
      \item \code{gene_name}: Character, gene name
      \item \code{CHR}: Numeric, chromosome
      \item \code{annotation}: Character, the Burden weight
      \item \code{category}: Character, the coding functional category used in the results of STAARpipeline. Choices include \code{plof}, \code{plof_ds}, \code{missense}, \code{disruptive_missense}, \code{synonymous}
      \item \code{Est}: Numeric, the Burden effect size estimates of gene-centric coding rare variant aggregate sets
      \item \code{Est_se}: Numeric, the standard error of Burden effect size estimates of gene-centric coding rare variant aggregate sets
      \item \code{pvalue}: Numeric, the P value of the Burden effect size estimates of gene-centric coding rare variant aggregate sets
    }
  }
  \item{noncoding}{A data frame with columns:
    \itemize{
      \item \code{gene_name}: Character, gene name
      \item \code{CHR}: Numeric, chromosome
      \item \code{annotation}: Character, the Burden weight
      \item \code{category}: Character, the coding functional category used in the results of STAARpipeline. Choices include \code{downstream}, \code{upstream}, \code{UTR}, \code{promoter_CAGE}, \code{promoter_DHS}, \code{enhancer_CAGE}, \code{enhancer_DHS}
      \item \code{Est}: Numeric, the Burden effect size estimates of gene-centric coding rare variant aggregate sets
      \item \code{Est_se}: Numeric, the standard error of Burden effect size estimates of gene-centric coding rare variant aggregate sets
      \item \code{pvalue}: Numeric, the P value of the Burden effect size estimates of gene-centric coding rare variant aggregate sets
    }
  }
  \item{window}{A data frame with columns:
    \itemize{
      \item \code{CHR}: Numeric, chromosome
      \item \code{annotation}: Character, the Burden weight
      \item \code{Start.Loc}: Character, the start location of sliding/dynamic window
      \item \code{End.Loc}: Character, the end location of sliding/dynamic window
      \item \code{Est}: Numeric, the Burden effect size estimates of gene-centric coding rare variant aggregate sets
      \item \code{Est_se}: Numeric, the standard error of Burden effect size estimates of gene-centric coding rare variant aggregate sets
      \item \code{pvalue}: Numeric, the P value of the Burden effect size estimates of gene-centric coding rare variant aggregate sets
    }
  }
}
}
\source{
The data is provided in \code{inst/extdata} of the \code{mr.carv} package, and loaded as a list containing two dataframe.
}
\usage{
TC
}
\description{
This dataset is stored in an Excel file \code{TC.xlsx} which contains two sheets:
\itemize{
  \item \code{indv_effect}: Individual variant associations with total cholesterol
  \item \code{coding}: Gene-Centric coding associations with total cholesterol
  \item \code{noncoding}: Gene-Centric noncoding associations with total cholesterol
  \item \code{window}: sliding/dynamic window region associations with total cholesterol
}
}
\examples{
# check the sample data from the package
str(TC$indv_effect)
str(TC$gene_coding)
}
\keyword{datasets}
