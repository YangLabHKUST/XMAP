comple <- function(allele) {
  ifelse(allele == "A", "T", ifelse(allele == "T", "A", ifelse(allele == "G", "C", ifelse(allele == "C", "G", allele))))
}

log_det <- function(X) {
  return(determinant(X)$modulus[1])
}
