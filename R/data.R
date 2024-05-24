#' TC Data
#'
#' This dataset is stored in an Excel file \code{TC.xlsx} which contains two sheets:
#' \itemize{
#'   \item \code{indv_effect}: Individual variant associations with total cholesterol
#'   \item \code{gene_effect}: Gene-Centric coding associations with total cholesterol
#' }
#'
#' @format A list with two data frames:
#' \describe{
#'   \item{indv_effect}{A data frame with columns:
#'     \itemize{
#'       \item \code{CHR}: Numeric, chromosome
#'       \item \code{POS}: Numeric, the position of each individual variant
#'       \item \code{REF}: Character, reference allele of each individual variant
#'       \item \code{ALT}: Character, alternative allele of each individual variant
#'       \item \code{annotation.id}: Character, the annotation id of each individual variant
#'       \item \code{Est}: Numeric, the effect size estimates of each individual variant on the exposure total cholesterol
#'       \item \code{Est_se}: Numeric, the standard error of the effect size estimates for each individual variant on the exposure total cholesterol
#'       \item \code{pvalue}: Numeric, the P value of the effect size for each individual variant
#'       \item \code{MAF}: Numeric, optional, for the harmonization
#'     }
#'   }
#'   \item{gene_effect}{A data frame with columns:
#'     \itemize{
#'       \item \code{gene_name}: Character, gene name
#'       \item \code{CHR}: Numeric, chromosome
#'       \item \code{annotation}: Character, the Burden weight
#'       \item \code{category}: Character, the coding functional category used in the results of STAARpipeline. Choices include \code{plof}, \code{plof_ds}, \code{missense}, \code{disruptive_missense}, \code{synonymous}
#'       \item \code{Est}: Numeric, the Burden effect size estimates of gene-centric coding rare variant aggregate sets
#'       \item \code{Est_se}: Numeric, the standard error of Burden effect size estimates of gene-centric coding rare variant aggregate sets
#'       \item \code{pvalue}: Numeric, the P value of the Burden effect size estimates of gene-centric coding rare variant aggregate sets
#'     }
#'   }
#' }
#' @source The data is provided in \code{inst/extdata} of the \code{mr.carv} package, and loaded as a list containing two dataframe.
#' @examples
#' # check the sample data from the package
#' str(TC$indv_effect)
#' str(TC$gene_effect)
#' @export
"TC"


# Ensure readxl package is available
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required to load the data.")
}

# Load the 'indv_effect' sheet
indv_effect <- readxl::read_excel(system.file("extdata", "TC.xlsx", package = "mr.carv"), sheet = "indv_effect")

# Load the 'gene_effect' sheet
gene_effect <- readxl::read_excel(system.file("extdata", "TC.xlsx", package = "mr.carv"), sheet = "gene_effect")

# Create a list to hold both sheets
TC <- list(
  indv_effect = indv_effect,
  gene_effect = gene_effect
)

