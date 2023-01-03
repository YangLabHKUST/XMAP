comple <- function(allele) {
  ifelse(allele == "A", "T", ifelse(allele == "T", "A", ifelse(allele == "G", "C", ifelse(allele == "C", "G", allele))))
}

log_det <- function(X) {
  return(determinant(X)$modulus[1])
}

get_pip <- function(X) {
  1 - apply(1 - X, 2, prod)
}


get_CS_XP <- function(object, Xcorr, component = c("shared", "pop1", "pop2"), squared = F, coverage = 0.95, min_abs_corr = 0.5, n_purity = 100) {
  
  gamma <- switch(component, "shared" = object$gamma, "pop1" = object$tgamma1, "pop2" = object$tgamma2)
  null_index = 0
  include_idx = rep(TRUE, nrow(gamma))
  
  # L x P binary matrix.
  status = in_CS(gamma, coverage)
  
  # L-list of CS positions.
  cs = lapply(1:nrow(status), function(i) which(status[i,] != 0))
  claimed_coverage = sapply(1:length(cs),
                            function(i) sum(gamma[i,][cs[[i]]]))
  
  include_idx = include_idx * (lapply(cs, length) > 0)
  include_idx = as.logical(include_idx)
  if (sum(include_idx) == 0)
    return(list(cs = NULL,
                coverage = NULL,
                requested_coverage = coverage))
  cs = cs[include_idx]
  claimed_coverage = claimed_coverage[include_idx]
  
  # prune by purity
  use_rfast = requireNamespace("Rfast", quietly = TRUE)
  purity = NULL
  for (i in 1:length(cs)) {
    if (null_index > 0 && null_index %in% cs[[i]])
      purity = rbind(purity, c(-9, -9, -9))
    else
      purity =
        rbind(purity,
              matrix(get_purity(cs[[i]], Xcorr, squared, n_purity, use_rfast), 1, 3))
  }
  purity = as.data.frame(purity)
  if (squared)
    colnames(purity) = c("min.sq.corr", "mean.sq.corr", "median.sq.corr")
  else
    colnames(purity) = c("min.abs.corr", "mean.abs.corr", "median.abs.corr")
  threshold = ifelse(squared, min_abs_corr^2, min_abs_corr)
  is_pure = which(purity[, 1] >= threshold)
  if (length(is_pure) > 0) {
    cs = cs[is_pure]
    purity = purity[is_pure,]
    row_names = paste0("L", which(include_idx)[is_pure])
    names(cs) = row_names
    rownames(purity) = row_names
    
    # Re-order CS list and purity rows based on purity.
    ordering = order(purity[, 1], decreasing = TRUE)
    return(list(cs = cs[ordering],
                purity = purity[ordering,],
                cs_index = which(include_idx)[is_pure[ordering]],
                coverage = claimed_coverage[ordering],
                requested_coverage = coverage))
  } else
    return(list(cs = NULL, coverage = NULL, requested_coverage = coverage))
}


get_CS <- function(object, X = NULL, Xcorr = NULL, squared = F, coverage = 0.95, min_abs_corr = 0.5, n_purity = 100, target_idx = 1:ncol(object$gamma)) {
  
  if(!is.null(X)) X <- X[,target_idx]
  if(!is.null(Xcorr)) Xcorr <- Xcorr[target_idx,target_idx]
  null_index = 0
  include_idx = rep(TRUE, nrow(object$gamma))
  
  # L x P binary matrix.
  status = in_CS(object, coverage, target_idx)
  
  # L-list of CS positions.
  cs = lapply(1:nrow(status), function(i) which(status[i,] != 0))
  claimed_coverage = sapply(1:length(cs),
                            function(i) sum(object$gamma[i,target_idx][cs[[i]]]))
  
  include_idx = include_idx * (lapply(cs, length) > 0)
  include_idx = as.logical(include_idx)
  if (sum(include_idx) == 0)
    return(list(cs = NULL,
                coverage = NULL,
                requested_coverage = coverage))
  cs = cs[include_idx]
  claimed_coverage = claimed_coverage[include_idx]
  
  # prune by purity
  use_rfast = requireNamespace("Rfast", quietly = TRUE)
  purity = NULL
  for (i in 1:length(cs)) {
    if (null_index > 0 && null_index %in% cs[[i]])
      purity = rbind(purity, c(-9, -9, -9))
    else
      purity =
        rbind(purity,
              matrix(get_purity(cs[[i]], X, Xcorr, squared, n_purity, use_rfast), 1, 3))
  }
  purity = as.data.frame(purity)
  if (squared)
    colnames(purity) = c("min.sq.corr", "mean.sq.corr", "median.sq.corr")
  else
    colnames(purity) = c("min.abs.corr", "mean.abs.corr", "median.abs.corr")
  threshold = ifelse(squared, min_abs_corr^2, min_abs_corr)
  is_pure = which(purity[, 1] >= threshold)
  if (length(is_pure) > 0) {
    cs = cs[is_pure]
    purity = purity[is_pure,]
    row_names = paste0("L", which(include_idx)[is_pure])
    names(cs) = row_names
    rownames(purity) = row_names
    
    # Re-order CS list and purity rows based on purity.
    ordering = order(purity[, 1], decreasing = TRUE)
    return(list(cs = cs[ordering],
                purity = purity[ordering,],
                cs_index = which(include_idx)[is_pure[ordering]],
                coverage = claimed_coverage[ordering],
                requested_coverage = coverage))
  } else
    return(list(cs = NULL, coverage = NULL, requested_coverage = coverage))
}

get_purity = function(pos, X, Xcorr, squared = FALSE, n = 100,
                      use_rfast) {
  if (missing(use_rfast))
    use_rfast = requireNamespace("Rfast", quietly = TRUE)
  if (use_rfast) {
    get_upper_tri = Rfast::upper_tri
    get_median = Rfast::med
  } else {
    get_upper_tri = function(R) R[upper.tri(R)]
    get_median = stats::median
  }
  if (length(pos) == 1)
    return(c(1, 1, 1))
  else {
    
    # Subsample the columns if necessary.
    if (length(pos) > n)
      pos = sample(pos, n)
    
    if (is.null(Xcorr)) {
      X_sub = X[, pos]
      X_sub = as.matrix(X_sub)
      value = abs(get_upper_tri(cor(X_sub)))
    } else
      value = abs(get_upper_tri(Xcorr[pos, pos]))
    if (squared)
      value = value^2
    if (squared)
      value = value^2
    return(c(min(value),
             sum(value) / length(value),
             get_median(value)))
  }
}


n_in_CS_x = function(x, coverage = 0.9)
  sum(cumsum(sort(x, decreasing = TRUE)) < coverage) + 1

# Return binary vector indicating if each point is in CS.
# x is a probability vector.
in_CS_x = function(x, coverage = 0.9) {
  n = n_in_CS_x(x, coverage)
  o = order(x, decreasing = TRUE)
  result = rep(0, length(x))
  result[o[1:n]] = 1
  return(result)
}

# Returns an l-by-p binary matrix indicating which variables are in
# susie credible sets.
in_CS = function(res, coverage = 0.9, target_idx) {
  
  res = res$gamma[,target_idx]
  return(t(apply(res, 1, function(x) in_CS_x(x, coverage))))
}


plot_CS <- function(pip, cs, b = NULL,...) {
  color = c(
    "dodgerblue2",
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )
  plot(pip, pch = 16, ylab = 'PIP', xlab = "variable",...)
  
  for (i in 1:length(cs)) {
    points(y = pip[cs[[i]]], x = cs[[i]], col = color[i], cex = 1.5, lwd = 2.5)
  }
  if (!is.null(b)) {
    points(which(b != 0), pip[b != 0], col = 2, pch = 16)
  }
}