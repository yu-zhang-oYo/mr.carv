# R/data.R

# Ensure readxl package is available
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required to load the data.")
}

# Load the 'indv_effect' sheet
indv_effect <- readxl::read_excel(system.file("extdata", "TC_CHR19.xlsx", package = "mr.carv"), sheet = "indv_effect")

# Load the 'gene_effect' sheet
gene_effect <- readxl::read_excel(system.file("extdata", "TC_CHR19.xlsx", package = "mr.carv"), sheet = "gene_effect")

# Create a list to hold both sheets
TC_CHR19 <- list(
  indv_effect = indv_effect,
  gene_effect = gene_effect
)

# Export the data
#' TC_CHR19 Data
#'
#' This dataset is stored in an Excel file `TC_CHR19.xlsx` which contains two sheets:
#' - `indv_effect`: Individual variant associations with total cholesterol
#' - `gene_effect`: Gene-Centric coding associations with total cholesterol
#'
#' @format A list with two data frames:
#' \describe{
#'   \item{indv_effect}{A data frame with columns:
#'     \itemize{
#'       \item{CHR}{Numeric, chromosome}
#'       \item{POS}{Numeric, the position of each individual variant}
#'       \item{REF}{Character, reference allele of each individual variant}
#'       \item{ALT}{Character, alternative allele of each individual variant}
#'       \item{annotation.id}{Character, the annotation id of each individual variant}
#'       \item{Est}{Numeric, the effect size estimates of each individual variant on the exposure total cholesterol}
#'       \item{Est_se}{Numeric, the standard error of the effect size estimates for each individual variant on the exposure total cholesterol}
#'       \item{pvalue}{Numeric, the P value of the effect size for each individual variant}
#'       \item{MAF}{Numeric, optional, for the harmonization}
#'     }
#'   }
#'   \item{gene_effect}{A data frame with columns:
#'     \itemize{
#'       \item{gene_name}{Character, gene name}
#'       \item{CHR}{Numeric, chromosome}
#'       \item{annotation}{Character, the Burden weight}
#'       \item{category}{Character, the coding functional category used in the results of STAARpipeline}
#'       \item{Est}{Numeric, the Burden effect size estimates of gene-centric coding rare variant aggregate sets}
#'       \item{Est_se}{Numeric, the standard error of Burden effect size estimates of gene-centric coding rare variant aggregate sets}
#'       \item{pvalue}{Numeric, the P value of the Burden effect size estimates of gene-centric coding rare variant aggregate sets}
#'     }
#'   }
#' }
#'
#' @source The data is provided as part of the `mr.carv` package.
#' @examples
#' # Load the sample data from the package
#' data(TC_CHR19)
#' str(TC_CHR19$indv_effect)
#' str(TC_CHR19$gene_effect)
#' @export
"TC_CHR19"
