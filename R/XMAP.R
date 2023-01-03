#'@title  Main function for XMAP
#'@description XMAP: Cross-population fine-mapping by leveraging genetic diversity and accounting for confounding bias
#'
#' @param R p by p by t array of reference LD matrices. p is the number of SNPs in the locus and t is the number of populations.
#' @param z p by t matrix of Z-scores from GWAS summary statistics.
#' @param n Vector of length t collecting the sample sizes of GWASs in the t populations.
#' @param K Maximum number of causal signals in the locus.
#' @param Sigma t by t by K array of the initial variance-covariance matrices of the K causal genetic effects.
#' @param Omega t by t variance-covariance matrix of the polygenic genetic effects. It should be estimated with the per-SNP heritability and per-SNP co-heritability by applying bi-variate LD score regression to genome-wide GWAS summary data.
#' @param Sig_E Vector of length t collecting the inflation constants of GWAS summary statistics in the t populations. It should be estimated with the intercepts of LD score regressions.
#' @param prior_weights Vector of length p collecting the final estimates of \code{pi}. Default is 1/p for all SNPs.
#' @param maxIter Maximum number of iterations in the VEM algorithm.
#' @param tol A non-negative number specifying the convergence tolerance of the relative difference between consecutive VEM iterations.
#' @param include_background A logical value indicating whether the polygenic  genetic effects should be included. If \code{Omega=NULL}, \code{include_background=FALSE}.
#' @param estimate_prior_variance A logical value indicating whether the prior variance \code{Sigma} should be updated in the VEM algorithm. Default is TRUE.
#' @param estimate_residual_variance A logical value indicating whether the prior variance \code{Sig_E} should be updated in the VEM algorithm. Default is FALSE. Note that this value should be always set to FALSE for identifiability.
#' @param estimate_background_variance A logical value indicating whether the prior variance \code{Omega} should be updated in the VEM algorithm. Default is FALSE. Note that this value should be always set to FALSE for identifiability.
#'
#'
#' @return a list with the following elements:
#' \describe{
#' \item{gamma: }{K by p matrix of posterior inclusion probability of the p SNPs in K causal signals.}
#' \item{mu: }{t by K by p array of posterior mean of causal effects.}
#' \item{nu: }{t by p matrix of posterior mean of polygenic effects.}
#' \item{Sigma: }{t by t by K array of variance-covariance matrices of the K causal genetic effects estimated by XMAP.}
#' \item{Omega: }{t by t variance-covariance matrix of the polygenic genetic effects.}
#' \item{Sig_E: }{Vector of length t collecting the inflation constants of GWAS summary statistics in the t populations.}
#' \item{ELBO: }{Vector of ELBO calues recorded at each step of VEM algorithm.}
#' \item{K: }{Number of causal signals.}
#' \item{prior_weights: }{Vector of length p collecting the final estimates of \code{pi}.}
#' }
#'
#' @examples
#' See vignette.
#' @export


XMAP <- function(R, z, n, K, Sigma = NULL, Omega = NULL, Sig_E = NULL, prior_weights = NULL, maxIter = 1000, tol = 1e-6, 
                 include_background = T, estimate_prior_variance = T, estimate_residual_variance = F, estimate_background_variance = F) {
  
  t <- dim(R)[3]
  p <- dim(R)[2]
  XX <- R
  for(tt in 1:t){
    XX[,,tt] <- XX[,,tt]*n[tt]
  }
  # XX <- aperm(apply(R, c(1,2), function(x) x*n),c(2,3,1))
  Xy <- t(t(z) * sqrt(n))
  if(is.null(Sigma)){
    Sigma <- array(replicate(K,diag(0.1,ncol = t,nrow = t)),dim = c(t,t,K))
    if(K>1){
      Sigma <- aperm(apply(Sigma, c(1,2), function(x) x/1:K),c(2,3,1))
    }
  }
  if(is.null(Omega)){
    include_background <- F
    cat("Background component Omega is not specified. Run XMAP-.\n")
  }
  if(is.null(Sig_E)){
    Sig_E <- rep(0.9,t)
  }
  if(is.null(prior_weights)){
    prior_weights <- rep(1/p,p)
  }
  
  if(include_background){
    fit_XMAP <- XMAP_suff(XX = XX,Xy = Xy,yty = n,n = n,K = K,Sigma = Sigma,Omega = Omega,Sig_E = Sig_E,prior_weights = prior_weights,
                          maxIter = maxIter,tol = tol,
                          estimate_prior_variance = estimate_prior_variance,estimate_residual_variance = estimate_residual_variance,
                          estimate_background_variance = estimate_background_variance)
  }

  return(fit_XMAP)
}
