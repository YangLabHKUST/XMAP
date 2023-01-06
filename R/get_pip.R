#' @title  Obtain the PIP of each SNP from XMAP output
#' @description This function computes the marginal PIP of each SNP from the gamma matrix of XMAP output
#'
#' @param X The gamma matrix from the XMAP output, a K by p matrix representing the PIP of p SNPs to be causal in K signals
#'
#' @return a vector of length p representing the PIP for each SNP to be causal.
#'
#' @examples
#' See vignette.
#' @export
get_pip <- function(X) {
  1 - apply(1 - X, 2, prod)
}




