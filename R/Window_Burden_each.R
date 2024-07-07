#' Get the Burden variable and Burden effect size for each sliding/dynamic window region
#'
#' The \code{Window_Burden_each} function takes in chromosome, start location, end location,
#' the object of opened annotated GDS file, and the object from fitting the null model to get the Burden variables and Burden Effect Size.
#' @param chr chromosome.
#' @param start_loc the start location of the sliding/dynamic window region that is significant in the existing summary statistics.
#' @param end_loc the end location of the sliding/dynamic window region that is significant in the existing summary statistics.
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function,
#' or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{\link{genesis2staar_nullmodel}} function.
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param rv_num_cutoff the cutoff of minimum number of variants of analyzing
#' a given variant-set (default = 2).
#' @param QC_label channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @param variant_type type of variant included in the analysis. Choices include "SNV", "Indel", or "variant" (default = "SNV").
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param Annotation_dir channel name of the annotations in the aGDS file \cr (default = "annotation/info/FunctionalAnnotation").
#' @param Annotation_name_catalog a data frame containing the name and the corresponding channel name in the aGDS file.
#' @param Beta_Annotation annotation weight used in the summary statistics.
#' @param silent logical: should the report of error messages be suppressed (default = FALSE).
#' @return a list containing the information of the given gene, Burden variable, and Burden effect size.
#' @export


Window_Burden_each <- function(chr,start_loc,end_loc,genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                  QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                  Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                  Beta_Annotation="Burden(1,1)",silent=FALSE){

  ## evaluate choices
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  phenotype.id <- as.character(obj_nullmodel$id_include)

  ## get SNV id
  filter <- seqGetData(genofile, QC_label)
  if(variant_type=="variant")
  {
    SNVlist <- filter == "PASS"
  }

  if(variant_type=="SNV")
  {
    SNVlist <- (filter == "PASS") & isSNV(genofile)
  }

  if(variant_type=="Indel")
  {
    SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  }

  variant.id <- seqGetData(genofile, "variant.id")

  ## Position
  position <- as.numeric(seqGetData(genofile, "position"))

  is.in <- (SNVlist)&(position>=start_loc)&(position<=end_loc)
  seqSetFilter(genofile,variant.id=variant.id[is.in],sample.id=phenotype.id)

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]

  # get the IDs
  annotation.id <- seqGetData(genofile, "annotation/id")
  variant.id <- seqGetData(genofile, "variant.id")

  ## impute missing
  if(!is.null(dim(Geno)))
  {
    if(dim(Geno)[2]>0)
    {
      if(geno_missing_imputation=="mean")
      {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor")
      {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }

  ## Annotation
  Annotation.PHRED <- NULL

  # Extract the parameters for the Burden function
  beta_par <- gsub("Burden\\(([^)]+)\\).*", "\\1", Beta_Annotation)
  # Split the parameters into separate values for beta function
  beta_par <- as.numeric(strsplit(beta_par, ",")[[1]])

  Annotation_name <- sub(".*?-(.*)", "\\1", Beta_Annotation)
  # Additional check if the extracted part is empty or same as input (meaning no actual extraction occurred)
  if(Annotation_name == Beta_Annotation) {
    Annotation_name <- NULL  # Set to NULL indicating no valid annotation after hyphen
  }

  if(variant_type=="SNV"){
    if(!is.null(Annotation_name)){
      # Check if the provided annotation name exists in the catalog
      if (!(Annotation_name %in% Annotation_name_catalog$name) & (Annotation_name != "aPC.LocalDiversity(-)")) {
        stop("The specified annotation name does not exist in the catalog.")
      }
      if(Annotation_name == "aPC.LocalDiversity(-)"){
        Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="aPC.LocalDiversity")]))
        Annotation.PHRED <- -10*log10(1-10^(-Annotation.PHRED/10))
      } else {
        Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name)]))
      }
      if(Annotation_name == "CADD"){
        Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
      }
    }
  }

  burden_effects <- 0
  # calculate the Burden effect size with Burden_EffectSize
  try(burden_effects <- Burden_Effect(genotype=Geno, obj_nullmodel=obj_nullmodel, beta_par=beta_par, annotation_phred=Annotation.PHRED, rare_maf_cutoff=rare_maf_cutoff, rv_num_cutoff=rv_num_cutoff),silent=silent)

  if(class(burden_effects)=="list")
  {
    # change the original code to output list
    results_temp <- list()
    results_temp$Chr <- chr
    results_temp$Start_loc <- start_loc
    results_temp$End_loc <- end_loc
    results_temp$Annotation <- ifelse(is.null(Annotation_name), paste0("Burden(",beta_par[1],",",beta_par[2],")"), paste0("Burden(",beta_par[1],",",beta_par[2],")-",Annotation_name))
    results_temp$'#SNV' <- burden_effects$num_variant

    # add the two kinds of IDs to the results
    results_temp$annotation.id <- annotation.id[burden_effects$RV_label]
    results_temp$variantIDs <- variant.id[burden_effects$RV_label]

    # add Burden effect size to the results
    results_temp$Burden_Effect_Size <- burden_effects$Burden_Effect_Size
    results_temp$Burden_Variable <- burden_effects$Burden_Variable
    rownames(results_temp$Burden_Variable) <- seqGetData(genofile, "sample.id")
  } else{
    results_temp <- NULL
  }

  seqResetFilter(genofile)

  return(results_temp)
}

