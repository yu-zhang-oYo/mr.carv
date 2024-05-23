# refer to the effect size calculation in metaSTAAR https://github.com/xihaoli/MetaSTAAR/blob/main/R/Burden_Burden_Effect_Size_meta.R
# Note: only check the function BurdenEffectSize_O_SMMAT because we have the data for this, still need to check the other two functions.

#' Calculate Burden Effect Size
#'
#' This function calculates the burden effect size for a given genotype matrix, considering relatedness of individuals and annotations.
#'
#' @param genotype A matrix or dgCMatrix of genotypes.
#' @param obj_nullmodel A null model object containing information about the model fit and residuals.
#' @param beta_par A vector of length 2 specifying the parameters of the beta distribution used for weighting. Default is \code{c(1, 1)}.
#' @param annotation_phred A numeric vector of annotation Phred scores. Default is \code{NULL}.
#' @param rare_maf_cutoff The cutoff for maximum minor allele frequency to define rare variants. Default is \code{0.01}.
#' @param rv_num_cutoff The minimum number of variants required in the set for analysis. Default is \code{2}.
#'
#' @return A list containing:
#' \item{num_variant}{The number of rare variants in the set.}
#' \item{cMAC}{The cumulative minor allele count of the variants.}
#' \item{RV_label}{A vector indicating which variants are rare.}
#' \item{Burden_Effect_Size}{A vector containing the burden effect size, standard error, score statistic, score standard error, and p-value.}
#' \item{Burden_Variable}{The burden variable weighted by the annotation score for each individual.}
#' @import Matrix
#' @importFrom stats pchisq
#' @export

# final calculation of effect size
# or we add this into the STAAR.R code so that it can output all the burden effect size and the Burden variable for certatin annotation,
# we need to revise the code in STAAR_pipeline to output the effect sizes
Burden_Effect <- function(genotype, obj_nullmodel, beta_par=c(1,1), annotation_phred = NULL, rare_maf_cutoff = 0.01, rv_num_cutoff = 2) {

  if(class(genotype)[1] != "matrix" && !(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix")){
    stop("genotype is not a matrix!")
  }

  if(dim(genotype)[2] == 1){
    stop("Number of rare variant in the set is less than 2!")
  }

  # annotation_phred <- as.data.frame(annotation_phred)
  if((!is.null(annotation_phred)) && dim(genotype)[2] != length(annotation_phred)){
    stop("Dimensions don't match for genotype and annotation!")
  }

  if(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix"){
    genotype <- as.matrix(genotype)
  }

  genotype <- matrix_flip(genotype)
  MAF <- genotype$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  Geno_rare <- genotype$Geno[, RV_label]

  rm(genotype)
  gc()

  if(sum(RV_label) >= rv_num_cutoff){
    G <- as(Geno_rare, "dgCMatrix")
    MAF <- MAF[RV_label]
    rm(Geno_rare)
    gc()

    # Burden weights
    ## default beta(1,1)
    w_B <- dbeta(MAF, beta_par[1], beta_par[2])

    if(!is.null(annotation_phred)){
      annotation_phred <- annotation_phred[RV_label]
      ## beta * the annotation function score
      annotation_rank <- 1 - 10^(-annotation_phred / 10)
      w_B <- annotation_rank * w_B
    }

    # weight the rare variants on a gene set into a Burden variable by annotation score for each individual
    # how to consider the relatedness of individuals when calculating the correlation of gene-sets
    Geno_rare_weighted <- G %*% w_B

    if(obj_nullmodel$relatedness){
      if(!obj_nullmodel$sparse_kins){
        P <- obj_nullmodel$P
        residuals.phenotype <- obj_nullmodel$scaled.residuals

        Burden_Effect_Size <- BurdenEffectSize_O_SMMAT(G, P, residuals=residuals.phenotype, weights_B=w_B)

      } else {
        Sigma_i <- obj_nullmodel$Sigma_i
        Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
        cov <- obj_nullmodel$cov

        residuals.phenotype <- obj_nullmodel$scaled.residuals

        Burden_Effect_Size <- BurdenEffectSize_O_SMMAT_sparse(G, Sigma_i, Sigma_iX, cov, residuals=residuals.phenotype, weights_B=w_B)
      }
    } else {
      X <- model.matrix(obj_nullmodel)
      working <- obj_nullmodel$weights
      sigma <- sqrt(summary(obj_nullmodel)$dispersion)
      if(obj_nullmodel$family[1] == "binomial"){
        fam <- 1
      }else if(obj_nullmodel$family[1] == "gaussian"){
        fam <- 0
      }

      residuals.phenotype <- obj_nullmodel$y - obj_nullmodel$fitted.values

      Burden_Effect_Size <- BurdenEffectSize_O(G, X, working, sigma, fam, residuals=residuals.phenotype, weights_B=w_B)
    }

    names(Burden_Effect_Size) <- c("Burden_Est", "Burden_SE_Est", "Burden_Score_Stat", "Burden_SE_Score", "Burden_pvalue")

    return(list(
      num_variant = sum(RV_label),
      cMAC = sum(G),
      RV_label = RV_label,
      Burden_Effect_Size = Burden_Effect_Size,
      Burden_Variable = Geno_rare_weighted
    ))
  } else {
    stop(paste0("Number of rare variant in the set is less than ", rv_num_cutoff, "!"))
  }
}



# the effect size under 3 situations
BurdenEffectSize_O <- function(G, X, working, sigma, fam, residuals, weights_B) {
  if (fam == 0) {
    tX_G <- t(X) %*% G
    Cov <- t(G) %*% G - t(tX_G) %*% solve(t(X) %*% X) %*% tX_G
  } else {
    tX_G <- t(X) %*% diag(working) %*% G
    Cov <- t(diag(working) %*% G) %*% G - t(tX_G) %*% solve(t(X) %*% diag(working) %*% X) %*% tX_G
  }
  x <- t(residuals) %*% G
  sum0 <- as.numeric(x %*% weights_B)
  sumw <- as.numeric(t(weights_B) %*% Cov %*% weights_B * sigma^2)
  burden_results <- c(sum0/sumw,
                      1/sqrt(sumw),
                      sum0,
                      sqrt(sumw),
                      pchisq(sum0^2/sumw, 1, lower.tail=FALSE))
  return(burden_results)
}

BurdenEffectSize_O_SMMAT <- function(G, P, residuals, weights_B) {
  Cov <- t(as.matrix(P %*% G)) %*% G  # P is not a matrix
  # Cov <- t(P %*% G) %*% G
  x <- t(residuals) %*% G
  sum0 <- as.numeric(x %*% weights_B)
  sumw <- as.numeric(t(weights_B) %*% Cov %*% weights_B)
  burden_results <- c(sum0/sumw,
                      1/sqrt(sumw),
                      sum0,
                      sqrt(sumw),
                      pchisq(sum0^2/sumw, 1, lower.tail=FALSE))
  return(burden_results)
}

BurdenEffectSize_O_SMMAT_sparse <- function(G, Sigma_i, Sigma_iX, cov, residuals, weights_B) {
  # Calculate the covariance matrix
  tSigma_iX_G <- t(Sigma_iX) %*% G
  Cov <- t(Sigma_i %*% G) %*% G - t(tSigma_iX_G) %*% cov %*% tSigma_iX_G
  x <- t(residuals) %*% G
  sum0 <- as.numeric(x %*% weights_B)
  sumw <- as.numeric(t(weights_B) %*% Cov %*% weights_B)
  burden_results <- c(sum0/sumw,
                      1/sqrt(sumw),
                      sum0,
                      sqrt(sumw),
                      pchisq(sum0^2/sumw, 1, lower.tail=FALSE))
  return(burden_results)
}





