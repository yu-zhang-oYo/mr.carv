#' Get the Burden variable and Burden effect size for each gene centric noncoding region
#'
#' The \code{Gene_Centric_Noncoding_Burden_each} function takes in chromosome, gene name, functional category,
#' the object of opened annotated GDS file, and the object from fitting the null model to get the Burden variables and Burden Effect Size.
#' @param chr chromosome.
#' @param gene_name name of the gene that is significant in the existing summary statistics.
#' @param category the coding functional category of \code{gene_name} that is significant in STAAR procedure. Choices include
#' \code{downstream}, \code{upstream}, \code{UTR}, \code{promoter_CAGE}, \code{promoter_DHS}, \code{enhancer_CAGE}, \code{enhancer_DHS}.
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



Gene_Centric_Noncoding_Burden_each <- function(chr, gene_name, category=c("downstream","upstream","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"),
                                            genofile, obj_nullmodel, rare_maf_cutoff=0.01, rv_num_cutoff=2,
                                            QC_label="annotation/filter", variant_type=c("SNV","Indel","variant"), geno_missing_imputation=c("mean","minor"),
                                            Annotation_dir="annotation/info/FunctionalAnnotation", Annotation_name_catalog,
                                            Beta_Annotation="Burden(1,1)",silent=FALSE){

  genes <- genes_info[genes_info[,2]==chr,]

  ## evaluate choices
  category <- match.arg(category)
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  phenotype.id <- as.character(obj_nullmodel$id_include)

  ## get SNV id, position, REF, ALT (whole genome)
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

  position <- as.numeric(seqGetData(genofile, "position"))
  variant.id <- seqGetData(genofile, "variant.id")

  rm(filter)
  gc()

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

  GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
  genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
  varid <- seqGetData(genofile, "variant.id")

  burden_effects <- 0

  ########################################
  #   Downstream

  if(category == "downstream"){
    is.in <- (GENCODE.Category=="downstream")&(SNVlist)
    variant.id.downstream <- variant.id[is.in]

    seqSetFilter(genofile,variant.id=variant.id.downstream,sample.id=phenotype.id)

    rm(variant.id.downstream)
    gc()

    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
    variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))

    variant.id.SNV <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)

    rm(GENCODE.Info)
    gc()

    rm(variant_gene_num)
    gc()

    Gene <- as.character(unlist(GENCODE.Info.split))

    rm(GENCODE.Info.split)
    gc()

    seqResetFilter(genofile)

    ### Gene
    is.in <- which(Gene==gene_name)
    variant.is.in <- variant.id.SNV[is.in]

    seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

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

    # get the rsID
    annotation.id <- seqGetData(genofile, "annotation/id")
    variant.id.gene <- seqGetData(genofile, "variant.id")

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

    # calculate the Burden effect size with Burden_EffectSize
    try(burden_effects <- Burden_Effect(genotype=Geno, obj_nullmodel=obj_nullmodel, beta_par=beta_par, annotation_phred=Annotation.PHRED, rare_maf_cutoff=rare_maf_cutoff, rv_num_cutoff=rv_num_cutoff),silent=silent)

  }

  ########################################
  #   Upstream

  if(category == "upstream"){
    is.in <- (GENCODE.Category=="upstream")&(SNVlist)
    variant.id.upstream <- variant.id[is.in]

    seqSetFilter(genofile,variant.id=variant.id.upstream,sample.id=phenotype.id)

    rm(variant.id.upstream)
    gc()

    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
    variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))

    variant.id.SNV <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)

    rm(GENCODE.Info)
    gc()

    rm(variant_gene_num)
    gc()

    Gene <- as.character(unlist(GENCODE.Info.split))

    rm(GENCODE.Info.split)
    gc()

    seqResetFilter(genofile)

    ### Gene
    is.in <- which(Gene==gene_name)
    variant.is.in <- variant.id.SNV[is.in]

    seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

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

    # get the rsID
    annotation.id <- seqGetData(genofile, "annotation/id")
    variant.id.gene <- seqGetData(genofile, "variant.id")

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

    # calculate the Burden effect size with Burden_EffectSize
    try(burden_effects <- Burden_Effect(genotype=Geno, obj_nullmodel=obj_nullmodel, beta_par=beta_par, annotation_phred=Annotation.PHRED, rare_maf_cutoff=rare_maf_cutoff, rv_num_cutoff=rv_num_cutoff),silent=silent)

  }


  ########################################################
  #                UTR

  if(category == "UTR"){
    is.in <- ((GENCODE.Category=="UTR3")|(GENCODE.Category=="UTR5")|(GENCODE.Category=="UTR5;UTR3"))&(SNVlist)
    variant.id.UTR <- variant.id[is.in]

    rm(GENCODE.Category)
    gc()

    seqSetFilter(genofile,variant.id=variant.id.UTR,sample.id=phenotype.id)

    rm(variant.id.UTR)
    gc()

    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[(]")

    rm(GENCODE.Info)
    gc()

    Gene <- as.character(sapply(GENCODE.Info.split,function(z) z[1]))

    rm(GENCODE.Info.split)
    gc()

    variant.id.SNV <- seqGetData(genofile, "variant.id")

    seqResetFilter(genofile)

    ### Gene
    is.in <- which(Gene==gene_name)
    variant.is.in <- variant.id.SNV[is.in]

    seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

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

    # get the rsID
    annotation.id <- seqGetData(genofile, "annotation/id")
    variant.id.gene <- seqGetData(genofile, "variant.id")

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

    # calculate the Burden effect size with Burden_EffectSize
    try(burden_effects <- Burden_Effect(genotype=Geno, obj_nullmodel=obj_nullmodel, beta_par=beta_par, annotation_phred=Annotation.PHRED, rare_maf_cutoff=rare_maf_cutoff, rv_num_cutoff=rv_num_cutoff),silent=silent)

  }


  #############################################
  #   Promoter-CAGE

  if(category == "Promoter_CAGE"){
    ## Promoter
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    promGobj <- promoters(genes(txdb), upstream = 3000, downstream = 3000)

    # Subsetting Promoters that within +/-3kb of TSS and have CAGE signals
    CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
    CAGEBvt <- CAGEAnno!=""
    CAGEidx <- which(CAGEBvt,useNames=TRUE)
    seqSetFilter(genofile,variant.id=varid[CAGEidx])
    seqSetFilter(genofile,promGobj,intersect=TRUE)
    CAGEpromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    CAGEGene <- unlist(lapply(strsplit(CAGEpromgene,"\\(|\\,|;|-"),`[[`,1))
    ##obtain variants info
    CAGEvchr <- as.numeric(seqGetData(genofile,"chromosome"))
    CAGEvpos <- as.numeric(seqGetData(genofile,"position"))
    CAGEvref <- as.character(seqGetData(genofile,"$ref"))
    CAGEvalt <- as.character(seqGetData(genofile,"$alt"))
    dfPromCAGEVarGene <- data.frame(CAGEvchr,CAGEvpos,CAGEvref,CAGEvalt,CAGEGene)

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
    variant.id.SNV <- variant.id[SNVlist]

    dfPromCAGEVarGene.SNV <- dfPromCAGEVarGene[SNVlist,]
    dfPromCAGEVarGene.SNV$CAGEvpos <- as.character(dfPromCAGEVarGene.SNV$CAGEvpos)
    dfPromCAGEVarGene.SNV$CAGEvref <- as.character(dfPromCAGEVarGene.SNV$CAGEvref)
    dfPromCAGEVarGene.SNV$CAGEvalt <- as.character(dfPromCAGEVarGene.SNV$CAGEvalt)

    seqResetFilter(genofile)

    rm(dfPromCAGEVarGene)
    gc()

    ### Gene
    is.in <- which(dfPromCAGEVarGene.SNV[,5]==gene_name)
    variant.is.in <- variant.id.SNV[is.in]

    seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

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

    # get the rsID
    annotation.id <- seqGetData(genofile, "annotation/id")
    variant.id.gene <- seqGetData(genofile, "variant.id")

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

    # calculate the Burden effect size with Burden_EffectSize
    try(burden_effects <- Burden_Effect(genotype=Geno, obj_nullmodel=obj_nullmodel, beta_par=beta_par, annotation_phred=Annotation.PHRED, rare_maf_cutoff=rare_maf_cutoff, rv_num_cutoff=rv_num_cutoff),silent=silent)

  }

  ##################################################
  #       Promoter-DHS

  if(category == "Promoter_DHS"){
    # Subsetting Promoters that within +/-3kb of TSS and have rOCRs signals
    rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
    rOCRsBvt <- rOCRsAnno!=""
    rOCRsidx <- which(rOCRsBvt,useNames=TRUE)
    seqSetFilter(genofile,variant.id=varid[rOCRsidx])

    seqSetFilter(genofile,promGobj,intersect=TRUE)
    rOCRspromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    rOCRsGene <- unlist(lapply(strsplit(rOCRspromgene,"\\(|\\,|;|-"),`[[`,1))
    ## obtain variants info
    rOCRsvchr <- as.numeric(seqGetData(genofile,"chromosome"))
    rOCRsvpos <- as.numeric(seqGetData(genofile,"position"))
    rOCRsvref <- as.character(seqGetData(genofile,"$ref"))
    rOCRsvalt <- as.character(seqGetData(genofile,"$alt"))
    dfPromrOCRsVarGene <- data.frame(rOCRsvchr,rOCRsvpos,rOCRsvref,rOCRsvalt,rOCRsGene)

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
    variant.id.SNV <- variant.id[SNVlist]

    dfPromrOCRsVarGene.SNV <- dfPromrOCRsVarGene[SNVlist,]
    dfPromrOCRsVarGene.SNV$rOCRsvpos <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvpos)
    dfPromrOCRsVarGene.SNV$rOCRsvref <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvref)
    dfPromrOCRsVarGene.SNV$rOCRsvalt <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvalt)

    seqResetFilter(genofile)

    rm(dfPromrOCRsVarGene)
    gc()

    ### Gene
    is.in <- which(dfPromrOCRsVarGene.SNV[,5]==gene_name)
    variant.is.in <- variant.id.SNV[is.in]

    seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

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

    # get the rsID
    annotation.id <- seqGetData(genofile, "annotation/id")
    variant.id.gene <- seqGetData(genofile, "variant.id")

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

    # calculate the Burden effect size with Burden_EffectSize
    try(burden_effects <- Burden_Effect(genotype=Geno, obj_nullmodel=obj_nullmodel, beta_par=beta_par, annotation_phred=Annotation.PHRED, rare_maf_cutoff=rare_maf_cutoff, rv_num_cutoff=rv_num_cutoff),silent=silent)

  }


  ###########################################
  #        Enhancer-CAGE

  if(category == "enhancer_CAGE"){
    #Now extract the GeneHancer with CAGE Signal Overlay
    genehancer <- genehancerAnno!=""

    CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
    CAGE <- CAGEAnno!=""
    CAGEGeneHancervt <- CAGEAnno!=""&genehancerAnno!=""
    CAGEGeneHanceridx <- which(CAGEGeneHancervt,useNames=TRUE)
    seqSetFilter(genofile,variant.id=varid[CAGEGeneHanceridx])

    # variants that covered by whole GeneHancer without CAGE overlap.
    genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
    enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
    enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
    enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
    enhancervpos <- as.numeric(seqGetData(genofile,"position"))
    enhancervref <- as.character(seqGetData(genofile,"$ref"))
    enhancervalt <- as.character(seqGetData(genofile,"$alt"))
    dfHancerVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)

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
    variant.id.SNV <- variant.id[SNVlist]

    dfHancerVarGene.SNV <- dfHancerVarGene[SNVlist,]
    dfHancerVarGene.SNV$enhancervpos <- as.character(dfHancerVarGene.SNV$enhancervpos)
    dfHancerVarGene.SNV$enhancervref <- as.character(dfHancerVarGene.SNV$enhancervref)
    dfHancerVarGene.SNV$enhancervalt <- as.character(dfHancerVarGene.SNV$enhancervalt)

    seqResetFilter(genofile)

    rm(dfHancerVarGene)
    gc()

    ### Gene
    is.in <- which(dfHancerVarGene.SNV[,5]==gene_name)
    variant.is.in <- variant.id.SNV[is.in]

    seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

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

    # get the rsID
    annotation.id <- seqGetData(genofile, "annotation/id")
    variant.id.gene <- seqGetData(genofile, "variant.id")

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

    # calculate the Burden effect size with Burden_EffectSize
    try(burden_effects <- Burden_Effect(genotype=Geno, obj_nullmodel=obj_nullmodel, beta_par=beta_par, annotation_phred=Annotation.PHRED, rare_maf_cutoff=rare_maf_cutoff, rv_num_cutoff=rv_num_cutoff),silent=silent)

  }

  if(category == "enhancer_DHS"){
    rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
    rOCRs <- rOCRsAnno!=""
    rOCRsGeneHancervt <- rOCRsAnno!=""&genehancerAnno!=""
    rOCRsGeneHanceridx <- which(rOCRsGeneHancervt,useNames=TRUE)
    seqSetFilter(genofile,variant.id=varid[rOCRsGeneHanceridx])
    # variants that covered by whole GeneHancer without rOCRs overlap.

    genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
    enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
    enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
    enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
    enhancervpos <- as.numeric(seqGetData(genofile,"position"))
    enhancervref <- as.character(seqGetData(genofile,"$ref"))
    enhancervalt <- as.character(seqGetData(genofile,"$alt"))
    dfHancerVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)

    rm(varid)
    gc()

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
    variant.id.SNV <- variant.id[SNVlist]

    dfHancerVarGene.SNV <- dfHancerVarGene[SNVlist,]
    dfHancerVarGene.SNV$enhancervpos <- as.character(dfHancerVarGene.SNV$enhancervpos)
    dfHancerVarGene.SNV$enhancervref <- as.character(dfHancerVarGene.SNV$enhancervref)
    dfHancerVarGene.SNV$enhancervalt <- as.character(dfHancerVarGene.SNV$enhancervalt)

    seqResetFilter(genofile)

    rm(dfHancerVarGene)
    gc()

    ### Gene
    is.in <- which(dfHancerVarGene.SNV[,5]==gene_name)
    variant.is.in <- variant.id.SNV[is.in]

    seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

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

    # get the rsID
    annotation.id <- seqGetData(genofile, "annotation/id")
    variant.id.gene <- seqGetData(genofile, "variant.id")

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

    # calculate the Burden effect size with Burden_EffectSize
    try(burden_effects <- Burden_Effect(genotype=Geno, obj_nullmodel=obj_nullmodel, beta_par=beta_par, annotation_phred=Annotation.PHRED, rare_maf_cutoff=rare_maf_cutoff, rv_num_cutoff=rv_num_cutoff),silent=silent)

  }


  if(class(burden_effects)=="list")
  {
    # change the original code to output list
    results_temp <- list()
    results_temp$Gene_name <- gene_name
    results_temp$Chr <- chr
    results_temp$Category <- category
    results_temp$Annotation <- ifelse(is.null(Annotation_name), paste0("Burden(",beta_par[1],",",beta_par[2],")"), paste0("Burden(",beta_par[1],",",beta_par[2],")-",Annotation_name))
    results_temp$'#SNV' <- burden_effects$num_variant

    # add the two kinds of IDs to the results
    results_temp$annotation.id <- annotation.id[burden_effects$RV_label]
    results_temp$variantIDs <- variant.id.gene[burden_effects$RV_label]

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
