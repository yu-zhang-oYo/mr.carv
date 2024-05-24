#' Estimate Individual Variant Effects
#'
#' This function estimates the effects of individual variants on the outcome, given the variant information, genotype data, and a null model.
#'
#' @param chr Chromosome number.
#' @param df_indv A data frame containing individual variant information with columns: `CHR`, `POS`, `annotation.id`, `REF`, and `ALT`.
#' @param genofile An opened annotated GDS (aGDS) file object containing genotype data.
#' @param obj_nullmodel A null model object containing information about the model fit and residuals.
#' @param mac_cutoff Minor allele count cutoff for rare variants. Default is 20.
#' @param QC_label Quality control label. Default is `"annotation/filter"`.
#' @param variant_type Type of variant to include in the analysis. Choices are `"variant"`, `"SNV"`, or `"Indel"`. Default is `"variant"`.
#' @param geno_missing_imputation Method for handling missing genotypes. Choices are `"mean"` or `"minor"`. Default is `"mean"`.
#'
#' @return A list containing:
#' \item{Variant_Estimate}{A data frame with columns: `CHR`, `variant.id`, `annotation.id`, `POS`, `REF`, `ALT`, `ALT_AF`, `MAF`, `N`, `pvalue`, `pvalue_log10`, `Score`, `Score_se`, `Est`, `Est_se`.}
#' \item{Variant_dosg}{A matrix of genotype dosages for the selected variants.}
#'
#' @import SeqArray
#' @importFrom SeqVarTools isSNV
#' @importFrom dplyr inner_join left_join
#' @export



# input the variant information to calculate the individual outcome to get the score estimates
Individual_Estimate <- function(chr, df_indv, genofile, obj_nullmodel, mac_cutoff=20,
                                QC_label="annotation/filter", variant_type=c("variant","SNV","Indel"), geno_missing_imputation=c("mean","minor")){

  ## evaluate choices
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  ## Null Model
  phenotype.id <- as.character(obj_nullmodel$id_include)

  samplesize <- length(phenotype.id)

  ### dense GRM
  if(!obj_nullmodel$sparse_kins)
  {
    P <- obj_nullmodel$P
    residuals.phenotype <- obj_nullmodel$scaled.residuals
  }

  ### sparse GRM
  if(obj_nullmodel$sparse_kins)
  {
    Sigma_i <- obj_nullmodel$Sigma_i
    Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
    cov <- obj_nullmodel$cov

    residuals.phenotype <- obj_nullmodel$scaled.residuals
  }

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

  is.in <- (SNVlist)
  variant.id <- seqGetData(genofile, "variant.id")
  variant.id <- variant.id[is.in]

  seqSetFilter(genofile,variant.id=variant.id,sample.id=phenotype.id)

  CHR <- as.numeric(seqGetData(genofile, "chromosome"))
  POS <- as.numeric(seqGetData(genofile, "position"))
  REF <- as.character(seqGetData(genofile, "$ref"))
  ALT <- as.character(seqGetData(genofile, "$alt"))
  annotation.id <- seqGetData(genofile, "annotation/id")
  variant.id <- seqGetData(genofile, "variant.id")

  # ? how to decide the condition here?
  # is.in <- (SNVlist)&(chr == CHR)&(pos == position)&(rs.id == rsIDs)
  # variant.id <- seqGetData(genofile, "variant.id")
  # variant.id <- variant.id[is.in]

  df_all <- data.frame(CHR=CHR, POS=POS, REF=REF, ALT=ALT, annotation.id=annotation.id, variant.id=variant.id)
  # df_indv <- data.frame(CHR=chrs, position=POSs, annotation.id=annotation.ids)
  # ? even in this way, there is still repeat rows in df_all
  # ? do not consider the harmony of the Reference allele.
  df_merge <- dplyr::inner_join(df_indv, df_all, by = c("CHR", "POS", "annotation.id", "REF", "ALT"))

  if(nrow(df_merge) == 0) {
    seqResetFilter(genofile)
    return(NULL)
  }

  seqSetFilter(genofile,variant.id=df_merge$variant.id, sample.id=phenotype.id)

  ## to make sure the Geno matrix have the same sample as in the outcome data, just in case some individuals in the outcome does not have the genotype.
  id.genotype <- seqGetData(genofile,"sample.id")

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]

  if(geno_missing_imputation=="mean")
  {
    Geno <- matrix_flip_mean(Geno)
  }
  if(geno_missing_imputation=="minor")
  {
    Geno <- matrix_flip_minor(Geno)
  }

  MAF <- Geno$MAF
  ALT_AF <- 1 - Geno$AF

  CHR <- as.numeric(seqGetData(genofile, "chromosome"))
  POS <- as.numeric(seqGetData(genofile, "position"))
  REF <- as.character(seqGetData(genofile, "$ref"))
  ALT <- as.character(seqGetData(genofile, "$alt"))
  N <- rep(samplesize,length(CHR))
  variant.id <- seqGetData(genofile, "variant.id")
  annotation.id <- seqGetData(genofile, "annotation/id")

  if(!all(CHR==chr))
  {
    warning("chr does not match the chromosome of genofile (the opened aGDS)!")
  }

  results <- c()
  Variant_dosg <- c()

  if(sum(MAF>=0.05)>=1)
  {
    ## Common_variants
    Geno_common <- Geno$Geno[,MAF>=0.05]

    CHR_common <- CHR[MAF>=0.05]
    POS_common <- POS[MAF>=0.05]
    REF_common <- REF[MAF>=0.05]
    ALT_common <- ALT[MAF>=0.05]
    MAF_common <- MAF[MAF>=0.05]
    ALT_AF_common <- ALT_AF[MAF>=0.05]
    N_common <- N[MAF>=0.05]

    variant.id_common <- variant.id[MAF >= 0.05]
    annotation.id_common <- annotation.id[MAF >= 0.05]

    if(sum(MAF>=0.05)==1)
    {
      Geno_common <- as.matrix(Geno_common,ncol=1)
    }

    ## sparse GRM
    if(obj_nullmodel$sparse_kins)
    {
      Score_test <- Individual_Score_Test(Geno_common, Sigma_i, Sigma_iX, cov, residuals.phenotype)
    }

    ## dense GRM
    if(!obj_nullmodel$sparse_kins)
    {
      Score_test <- Individual_Score_Test_denseGRM(Geno_common, P, residuals.phenotype)
    }

    results_temp <- data.frame(CHR=CHR_common,variant.id=variant.id_common, annotation.id=annotation.id_common,
                               POS=POS_common,REF=REF_common,ALT=ALT_common,ALT_AF=ALT_AF_common,MAF=MAF_common,N=N_common,
                               pvalue=exp(-Score_test$pvalue_log),pvalue_log10=Score_test$pvalue_log/log(10),
                               Score=Score_test$Score,Score_se=Score_test$Score_se,
                               Est=Score_test$Est,Est_se=Score_test$Est_se)
    results <- rbind(results,results_temp)
    Variant_dosg <- cbind(Variant_dosg, Geno_common)
  }

  ## Rare_variants
  if(sum((MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05))>=1)
  {
    Geno_rare <- Geno$Geno[,(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
    Variant_dosg <- cbind(Variant_dosg, Geno_rare)
    CHR_rare <- CHR[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
    POS_rare <- POS[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
    REF_rare <- REF[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
    ALT_rare <- ALT[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
    MAF_rare <- MAF[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
    ALT_AF_rare <- ALT_AF[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
    N_rare <- N[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]

    # get the id number
    variant.id_rare <- variant.id[(MAF > (mac_cutoff - 0.5)/samplesize/2) & (MAF < 0.05)]
    annotation.id_rare <- annotation.id[(MAF > (mac_cutoff - 0.5)/samplesize/2) & (MAF < 0.05)]

    ## sparse GRM
    if(obj_nullmodel$sparse_kins)
    {
      if(sum((MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05))>=2)
      {
        Geno_rare <- as(Geno_rare,"dgCMatrix")
        Score_test <- Individual_Score_Test_sp(Geno_rare, Sigma_i, Sigma_iX, cov, residuals.phenotype)
      }else
      {
        Geno_rare <- as.matrix(Geno_rare,ncol=1)
        Score_test <- Individual_Score_Test(Geno_rare, Sigma_i, Sigma_iX, cov, residuals.phenotype)
      }
    }

    ## dense GRM
    if(!obj_nullmodel$sparse_kins)
    {
      if(sum((MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05))>=2)
      {
        Geno_rare <- as(Geno_rare,"dgCMatrix")
        Score_test <- Individual_Score_Test_sp_denseGRM(Geno_rare, P, residuals.phenotype)
      }else
      {
        Geno_rare <- as.matrix(Geno_rare,ncol=1)
        Score_test <- Individual_Score_Test_denseGRM(Geno_rare, P, residuals.phenotype)
      }
    }

    results_temp <- data.frame(CHR=CHR_rare, variant.id=variant.id_rare, annotation.id=annotation.id_rare,
                               POS=POS_rare,REF=REF_rare,ALT=ALT_rare,ALT_AF=ALT_AF_rare,MAF=MAF_rare,N=N_rare,
                               pvalue=exp(-Score_test$pvalue_log),pvalue_log10=Score_test$pvalue_log/log(10),
                               Score=Score_test$Score,Score_se=Score_test$Score_se,
                               Est=Score_test$Est,Est_se=Score_test$Est_se)

    results <- rbind(results,results_temp)
  }

  if(!is.null(Variant_dosg)){
    rownames(Variant_dosg) <- phenotype.id
  } else{
    seqResetFilter(genofile)
    return(NULL)
  }

  seqResetFilter(genofile)
  return(list(Variant_Estimate=results, Variant_dosg=Variant_dosg))
}


