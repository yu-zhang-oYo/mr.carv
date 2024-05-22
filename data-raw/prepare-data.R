# data-raw/prepare-data.R

library(readxl)
library(usethis)

# Read the 'indv_effect' and 'gene_effect' sheets from the Excel file
file_path <- "inst/extdata/TC_CHR19.xlsx"

indv_effect <- read_excel(file_path, sheet = "indv_effect")
gene_effect <- read_excel(file_path, sheet = "gene_effect")

# Combine the data if needed, otherwise save as separate objects
TC_CHR19 <- list(indv_effect = indv_effect, gene_effect = gene_effect)

# # Read the individual analysis data
# individual_analysis_path <- "inst/extdata/individual_analysis_19.out"
# individual_analysis_19 <- read.table(individual_analysis_path, header = TRUE, quote = "", stringsAsFactors = FALSE)

# Save the datasets in the data directory
usethis::use_data(TC_CHR19, overwrite = TRUE)
