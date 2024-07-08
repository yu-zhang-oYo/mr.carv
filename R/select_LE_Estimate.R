#' Select the genes and single variant that is in linkage equilibrium
#'
#' The \code{select_LE_Estimate} function finds the individual variants and gene-centric variants that are in linkage equilibrium
#' on a given chromosome. It calculates the correlation between genes and individual SNPs and selects a maximal independent set
#' of SNPs.
#'
#' @param chr Chromosome number.
#' @param df_indv Dataframe containing individual variant information with columns: CHR, POS, annotation.id, REF, ALT, Est, Est_se, and pvalue.
#' @param df_coding Dataframe containing coding gene information with columns: CHR, gene_name, category, annotation, Est, Est_se, and pvalue.
#' @param df_noncoding Dataframe containing noncoding gene information with columns: CHR, gene_name, category, annotation, Est, Est_se, and pvalue.
#' @param df_window Dataframe containing window region information with columns: CHR, Start.Loc, End.Loc, annotation, Est, Est_se, and pvalue.
#' @param genofile An object of an opened annotated GDS (aGDS) file.
#' @param obj_nullmodel An object from fitting the null model, which is either the output from the \code{\link{fit_nullmodel}} function,
#' or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{\link{genesis2staar_nullmodel}} function.
#' @param rare_maf_cutoff Minor allele frequency cutoff for rare variants. Default is 0.01.
#' @param rv_num_cutoff Cutoff for the minimum number of variants in a given variant-set. Default is 2.
#' @param QC_label Channel name of the QC label in the GDS/aGDS file. Default is "annotation/filter".
#' @param variant_type_indv Type of individual variants to include in the analysis. Choices include "SNV", "Indel", or "variant". Default is c("SNV", "Indel", "variant").
#' @param variant_type_gene Type of gene-centric variants to include in the analysis. Choices include "SNV", "Indel", or "variant". Default is c("SNV", "Indel", "variant").
#' @param geno_missing_imputation Method for imputing missing genotypes. Choices are "mean" or "minor". Default is c("mean", "minor").
#' @param Annotation_dir Directory for functional annotation information in the aGDS file. Default is "annotation/info/FunctionalAnnotation".
#' @param Annotation_name_catalog A data frame containing the name and the corresponding channel name in the aGDS file.
#' @param le_threshold Threshold for linkage disequilibrium (LD) correlation to define independent sets. Default is 0.1.
#' @param pvalues_weight Logical; if TRUE, uses p-values to weight the selection of independent SNPs, if FALSE, randomly select one set with independent SNPs. Default is FALSE.
#' @param silent Logical; if TRUE, suppresses messages. Default is FALSE.
#'
#' @return A list containing:
#' \describe{
#'   \item{X_indv}{Dataframe of individual variants with selected independent SNPs for exposure.}
#'   \item{X_gene}{Dataframe of gene-centric variants with selected independent SNPs for exposure.}
#'   \item{Y_indv}{Dataframe of effect size estimates for individual variants for outcome.}
#'   \item{Y_gene}{Dataframe of effect size estimates for gene-centric variants for outcome.}
#' }
#' @export

select_LE_Estimate <- function(chr, df_indv, df_coding, df_noncoding, df_window, genofile, obj_nullmodel, rare_maf_cutoff=0.01, rv_num_cutoff=2,
                               QC_label="annotation/filter",variant_type_indv=c("SNV","Indel","variant"),variant_type_gene=c("SNV","Indel","variant"),
                               geno_missing_imputation=c("mean","minor"),
                               Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                               le_threshold=0.1, pvalues_weight=FALSE, silent=FALSE){

  # Get the individual variant summary effect on the chromosome
  df_indv_chr <- df_indv[df_indv$CHR == chr, ]
  if (nrow(df_indv_chr) == 0) {
    results_indv <- NULL
  } else {
    results_indv <- Individual_Estimate(chr=chr, df_indv=df_indv_chr, genofile=genofile, obj_nullmodel=obj_nullmodel,
                                        QC_label=QC_label, variant_type=variant_type_indv, geno_missing_imputation=geno_missing_imputation)
  }

  if(!is.null(results_indv)){
    # Keep the rows in df_indv_chr that match results_indv[["Variant_Estimate"]]
    X_indv <- merge(df_indv_chr, results_indv[["Variant_Estimate"]][,c("CHR", "POS", "annotation.id", "REF", "ALT")], by = c("CHR", "POS", "annotation.id", "REF", "ALT"))
    X_indv <- X_indv[order(match(X_indv$annotation.id, results_indv[["Variant_Estimate"]]$annotation.id)), ]
  } else{
    X_indv <- NULL
  }

  # Process coding regions
  df_coding_chr <- df_coding[df_coding$CHR == chr, ]
  if (nrow(df_coding_chr) == 0) {
    results_coding <- NULL
  } else {
    results_coding <- Gene_Centric_Coding_Burden(chr=chr, df_gene=df_coding_chr, genofile=genofile, obj_nullmodel=obj_nullmodel,
                                                 rare_maf_cutoff=rare_maf_cutoff, QC_label=QC_label,
                                                 variant_type=variant_type_gene, geno_missing_imputation=geno_missing_imputation,
                                                 Annotation_dir=Annotation_dir, Annotation_name_catalog=Annotation_name_catalog, silent=silent)
  }

  if(!is.null(results_coding)){
    results_coding_chr <- lapply(results_coding, function(x) {
      data.frame(
        CHR = x$Chr,
        gene_name = x$Gene_name,
        category = x$Category,
        annotation = x$Annotation,
        Est = x$Burden_Effect_Size[1],
        Est_se = x$Burden_Effect_Size[2],
        pvalue = x$Burden_Effect_Size[5],
        "No.of SNV" = x$`#SNV`
      )
    }) %>% do.call(rbind, .)
    # Add a temporary ID to results_coding_chr to preserve the order of results
    results_coding_chr$temp_id <- 1:nrow(results_coding_chr)
    X_coding <- merge(df_coding_chr, results_coding_chr[,c("CHR", "gene_name", "category", "annotation", "temp_id")], by = c("CHR", "gene_name", "category", "annotation"))
    X_coding <- X_coding[order(match(X_coding$temp_id, results_coding_chr$temp_id)), ]
    X_coding$temp_id <- NULL
  } else {
    X_coding <- NULL
    results_coding_chr <- NULL
  }

  # Process noncoding regions
  df_noncoding_chr <- df_noncoding[df_noncoding$CHR == chr, ]
  if (nrow(df_noncoding_chr) == 0) {
    results_noncoding <- NULL
  } else {
    results_noncoding <- Gene_Centric_Noncoding_Burden(chr=chr, df_gene=df_noncoding_chr, genofile=genofile, obj_nullmodel=obj_nullmodel,
                                                       rare_maf_cutoff=rare_maf_cutoff, rv_num_cutoff=rv_num_cutoff,
                                                       QC_label=QC_label, variant_type=variant_type_gene, geno_missing_imputation=geno_missing_imputation,
                                                       Annotation_dir=Annotation_dir, Annotation_name_catalog=Annotation_name_catalog, silent=silent)
  }

  if (!is.null(results_noncoding)) {
    results_noncoding_chr <- lapply(results_noncoding, function(x) {
      data.frame(
        CHR = x$Chr,
        gene_name = x$Gene_name,
        category = x$Category,
        annotation = x$Annotation,
        Est = x$Burden_Effect_Size[1],
        Est_se = x$Burden_Effect_Size[2],
        pvalue = x$Burden_Effect_Size[5],
        "No.of SNV" = x$`#SNV`
      )
    }) %>% do.call(rbind, .)
    # Add a temporary ID to results_noncoding_chr to preserve the order of results
    results_noncoding_chr$temp_id <- 1:nrow(results_noncoding_chr)
    X_noncoding <- merge(df_noncoding_chr, results_noncoding_chr[, c("CHR", "gene_name", "category", "annotation", "temp_id")], by = c("CHR", "gene_name", "category", "annotation"))
    X_noncoding <- X_noncoding[order(match(X_noncoding$temp_id, results_noncoding_chr$temp_id)), ]
    X_noncoding$temp_id <- NULL
  } else {
    X_noncoding <- NULL
    results_noncoding_chr <- NULL
  }

  # Process window regions
  df_window_chr <- df_window[df_window$CHR == chr, ]
  if (nrow(df_window_chr) == 0) {
    results_window <- NULL
  } else {
    results_window <- Window_Burden(chr=chr, df_window=df_window_chr, genofile=genofile, obj_nullmodel=obj_nullmodel,
                                    rare_maf_cutoff=rare_maf_cutoff, rv_num_cutoff=rv_num_cutoff,
                                    QC_label=QC_label, variant_type=variant_type_gene, geno_missing_imputation=geno_missing_imputation,
                                    Annotation_dir=Annotation_dir, Annotation_name_catalog=Annotation_name_catalog, silent=silent)
  }

  if (!is.null(results_window)) {
    results_window_chr <- lapply(results_window, function(x) {
      data.frame(
        CHR = x$Chr,
        Start.Loc = x$Start_loc,
        End.Loc = x$End_loc,
        annotation = x$Annotation,
        Est = x$Burden_Effect_Size[1],
        Est_se = x$Burden_Effect_Size[2],
        pvalue = x$Burden_Effect_Size[5],
        "No.of SNV" = x$`#SNV`
      )
    }) %>% do.call(rbind, .)
    # Add a temporary ID to results_window_chr to preserve the order of results
    results_window_chr$temp_id <- 1:nrow(results_window_chr)
    X_window <- merge(df_window_chr, results_window_chr[, c("CHR", "Start.Loc", "End.Loc", "annotation", "temp_id")], by = c("CHR", "Start.Loc", "End.Loc", "annotation"))
    X_window <- X_window[order(match(X_window$temp_id, results_window_chr$temp_id)), ]
    X_window$temp_id <- NULL
  } else {
    X_window <- NULL
    results_window_chr <- NULL
  }

  # Calculate correlation between genes and individual SNPs
  Burden_Variables <- NULL
  if (!is.null(results_coding)) {
    Burden_Variables <- as.matrix(do.call(cbind, lapply(results_coding, function(x) x$Burden_Variable)))
  }
  if (!is.null(results_noncoding)) {
    Burden_Variables_noncoding <- as.matrix(do.call(cbind, lapply(results_noncoding, function(x) x$Burden_Variable)))
    Burden_Variables <- cbind(Burden_Variables, Burden_Variables_noncoding)
  }
  if (!is.null(results_window)) {
    Burden_Variables_window <- as.matrix(do.call(cbind, lapply(results_window, function(x) x$Burden_Variable)))
    Burden_Variables <- cbind(Burden_Variables, Burden_Variables_window)
  }
  if (!is.null(Burden_Variables)) {
    Burden_Variables <- as.data.frame(Burden_Variables)
  }

  if (!is.null(results_indv)) {
    dosg_indv <- as.data.frame(results_indv$Variant_dosg)
    if (!is.null(Burden_Variables)) {
      # Merge by sample id
      Burden_Variables <- merge(Burden_Variables, dosg_indv, by = "row.names", all = FALSE)
      Burden_Variables <- Burden_Variables[,-1]
    } else {
      Burden_Variables <- dosg_indv
    }
  }

  if (!is.null(Burden_Variables)) {
    colnames(Burden_Variables) <- 1:ncol(Burden_Variables)
    cor_Burden_Variables <- cor(Burden_Variables, use="pairwise.complete.obs", method = "pearson")

    # Convert the correlation matrix to an adjacency matrix for LD
    adj_matrix <- abs(cor_Burden_Variables) > le_threshold
    diag(adj_matrix) <- FALSE  # Remove self-loops

    # Create a graph from the adjacency matrix
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")

    # Find a maximal independent set
    independent_sets <- largest_ivs(g)
    if (!pvalues_weight) {
      selected_indx <- unlist(sample(independent_sets, 1))
    } else {
      # Get p-values from the original data frames
      pvalues_coding <- if (!is.null(X_coding)) X_coding$pvalue else NULL
      pvalues_noncoding <- if (!is.null(X_noncoding)) X_noncoding$pvalue else NULL
      pvalues_window <- if (!is.null(X_window)) X_window$pvalue else NULL
      pvalues_indv <- if (!is.null(X_indv)) X_indv$pvalue else NULL

      # Combine p-values from genes and individuals
      pvalues <- c(pvalues_coding, pvalues_noncoding, pvalues_window, pvalues_indv)

      # Using lapply to calculate the sum/mean of P-values in each set
      sums_pvalues <- lapply(independent_sets, function(x) {
        mean(pvalues[x])  # Sum the P-values at these indices
      })

      # Find the index of the set with the minimum sum of P-values
      min_index <- which.min(unlist(sums_pvalues))
      selected_indx <- independent_sets[[min_index]]
    }

    num_coding <- if (!is.null(results_coding)) length(results_coding) else 0
    num_noncoding <- if (!is.null(results_noncoding)) length(results_noncoding) else 0
    num_window <- if (!is.null(results_window)) length(results_window) else 0

    selected_coding_indx <- selected_indx[selected_indx <= num_coding]
    selected_noncoding_indx <- as.numeric(selected_indx[(selected_indx > num_coding) & (selected_indx <= num_coding + num_noncoding)]) - num_coding
    selected_window_indx <- as.numeric(selected_indx[(selected_indx > num_coding + num_noncoding) & (selected_indx <= num_coding + num_noncoding + num_window)]) - (num_coding + num_noncoding)
    selected_indv_indx <- as.numeric(selected_indx[selected_indx > num_coding + num_noncoding + num_window]) - (num_coding + num_noncoding + num_window)

    results_coding_chr <- if ((!is.null(results_coding_chr)) & (length(selected_coding_indx)>0)) results_coding_chr[selected_coding_indx, ] else NULL
    results_noncoding_chr <- if ((!is.null(results_noncoding_chr)) & (length(selected_noncoding_indx)>0)) results_noncoding_chr[selected_noncoding_indx, ] else NULL
    results_window_chr <- if ((!is.null(results_window_chr)) & (length(selected_window_indx)>0)) results_window_chr[selected_window_indx, ] else NULL
    results_indv <- if ((!is.null(results_indv)) & (length(selected_indv_indx)>0)) results_indv$Variant_Estimate[selected_indv_indx,] else NULL

    X_indv <- if(!is.null(results_indv)) X_indv[selected_indv_indx, ] else NULL
    X_coding <- if(!is.null(results_coding_chr)) X_coding[selected_coding_indx, ] else NULL
    X_noncoding <- if(!is.null(results_noncoding_chr)) X_noncoding[selected_noncoding_indx, ] else NULL
    X_window <- if(!is.null(results_window_chr)) X_window[selected_window_indx, ] else NULL

    return(list(X_indv=X_indv, X_coding=X_coding, X_noncoding=X_noncoding, X_window=X_window, Y_indv = results_indv, Y_coding = results_coding_chr, Y_noncoding = results_noncoding_chr, Y_window = results_window_chr))
  } else {
    return(list(X_indv=NULL, X_coding=NULL, X_noncoding=NULL, X_window=NULL, Y_indv=NULL, Y_gene=NULL, Y_noncoding=NULL, Y_window=NULL))
  }
}
