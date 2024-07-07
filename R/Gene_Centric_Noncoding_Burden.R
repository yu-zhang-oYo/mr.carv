#' Get the Burden variable and Burden effect size of the gene centric noncoding region
#'
#' The \code{Gene_Centric_Noncoding_Burden} function performs gene-centric noncoding burden analysis for each gene in the given dataframe.
#'
#' @param chr Chromosome number.
#' @param df_gene Dataframe containing gene information and genes' association result with the exposure. Columns include CHR, gene_name, category, annotation, Est, Est_se, pvalue.
#' @param genofile An object of an opened annotated GDS (aGDS) file.
#' @param obj_nullmodel An object from fitting the null model, which is either the output from the \code{\link{fit_nullmodel}} function,
#' or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{\link{genesis2staar_nullmodel}} function.
#' @param rare_maf_cutoff Minor allele frequency cutoff for rare variants. Default is 0.01.
#' @param rv_num_cutoff Cutoff for the minimum number of variants in a given variant-set. Default is 2.
#' @param QC_label Channel name of the QC label in the GDS/aGDS file. Default is "annotation/filter".
#' @param variant_type Type of variants to include in the analysis. Choices include "SNV", "Indel", or "variant". Default is c("SNV", "Indel", "variant").
#' @param geno_missing_imputation Method for imputing missing genotypes. Choices are "mean" or "minor". Default is c("mean", "minor").
#' @param Annotation_dir Directory for functional annotation information in the aGDS file. Default is "annotation/info/FunctionalAnnotation".
#' @param Annotation_name_catalog A data frame containing the name and the corresponding channel name in the aGDS file.
#' @param silent Logical; if TRUE, suppresses messages. Default is FALSE.
#'
#' @return A list of results from gene-centric noncoding burden analysis for each gene.
#' @export



Gene_Centric_Noncoding_Burden <- function(chr, df_gene, genofile, obj_nullmodel, rare_maf_cutoff=0.01, rv_num_cutoff=2,
                                       QC_label="annotation/filter", variant_type=c("SNV","Indel","variant"), geno_missing_imputation=c("mean","minor"),
                                       Annotation_dir="annotation/info/FunctionalAnnotation", Annotation_name_catalog, silent=FALSE){
  results <- list()
  for (i in 1:nrow(df_gene)) {
    results[[i]] <- Gene_Centric_Noncoding_Burden_each(chr=df_gene$CHR[i], gene_name=df_gene$gene_name[i], category=df_gene$category[i],
                                                    genofile=genofile, obj_nullmodel=obj_nullmodel, rare_maf_cutoff=rare_maf_cutoff, rv_num_cutoff=rv_num_cutoff,
                                                    QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                    Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                    Beta_Annotation= df_gene$annotation[i],silent=FALSE)
  }
  if(length(results)==0){
    results <- NULL
  } else {
    results <- Filter(Negate(is.null), results)
  }
  return(results)
}

