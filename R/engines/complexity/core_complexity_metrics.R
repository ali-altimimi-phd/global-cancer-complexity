#' This script defines the core structural metrics used to quantify transcriptomic
#' complexity, treating complexity as a function of matrix conditioning, spectral
#' organization, effective dimensionality, and related summaries.
#' 
#' In this framework, malignant transformation is evaluated in terms of how it 
#' alters the geometry and stability of the expression matrix rather than 
#' only mean-level expression changes.
NULL
# -------------------------------------------------------------
# Core Complexity Metric Functions
# -------------------------------------------------------------

#' Compute covariance matrix condition number
#'
#' @param mat A numeric matrix (samples in columns).
#' @return Condition number of the covariance matrix (rounded).
#' @export
get_cov_kappa <- function(mat) {
  mat <- mat[complete.cases(mat), ]
  cov_matrix <- cov(mat)
  round(pracma::cond(cov_matrix), 2)
}

#' Compute 2-norm condition number
#'
#' @param mat A numeric matrix.
#' @return 2-norm condition number (rounded).
#' @export
get_2norm_kappa <- function(mat) {
  mat <- mat[complete.cases(mat), ]
  round(kappa(mat), 2)
}

#' Compute SVD-based condition number
#'
#' @param mat A numeric matrix.
#' @return Ratio of largest to smallest singular values (rounded).
#' @export
get_svd_kappa <- function(mat) {
  if (!is.matrix(mat)) return(NA_real_)
  mat <- mat[complete.cases(mat), ]
  svd_d <- svd(mat)$d
  round(max(svd_d) / max(min(svd_d), .Machine$double.eps), 2)
}

#' Calculate effective rank based on entropy of eigenvalues
#'
#' @param mat A numeric matrix.
#' @return Effective rank as exp(entropy of eigenvalue distribution).
#' @export
get_effective_rank <- function(mat) {
  if (!is.matrix(mat)) return(NA_real_)
  mat <- mat[complete.cases(mat), ]
  eig_vals <- eigen(cov(mat), only.values = TRUE)$values
  eig_vals <- eig_vals[eig_vals > 0]
  p <- eig_vals / sum(eig_vals)
  h <- -sum(p * log(p))
  round(exp(h), 2)
}

#' Estimate matrix sparsity
#'
#' @param mat A numeric matrix.
#' @param threshold Threshold for considering values as zero.
#' @return Proportion of near-zero elements.
#' @export
get_matrix_sparsity <- function(mat, threshold = 1e-5) {
  if (!is.matrix(mat)) return(NA_real_)
  mat <- mat[complete.cases(mat), ]
  round(mean(abs(mat) < threshold), 3)
}

#' Compute composite kappa from three condition number estimates
#'
#' @param mat A numeric matrix.
#' @return Mean of three condition numbers.
#' @export
compute_composite_kappa <- function(mat) {
  if (!is.matrix(mat)) return(NA_real_)
  mat <- mat[complete.cases(mat), ]
  k1 <- get_cov_kappa(mat)
  k2 <- get_2norm_kappa(mat)
  k3 <- get_svd_kappa(mat)
  round(mean(c(k1, k2, k3), na.rm = TRUE), 2)
}

#' Determine complexity direction (gain/loss)
#'
#' @param diff Difference in complexity between tumor and normal.
#' @return "gained", "lost", or NA.
#' @export
complexity_gain <- function(diff) {
  if (is.na(diff)) return(NA_character_)
  if (diff > 0) "gained" else "lost"
}

#' Main complexity wrapper function
#'
#' @param tissue Tissue name or label.
#' @param mat A numeric matrix.
#' @return A list of complexity metrics.
#' @export
main_complexity <- function(tissue, mat) {
  if (is.null(mat) || !is.matrix(mat)) {
    return(list(
      tissue     = tissue, samples = NA,
      kappa1     = NA, kappa2 = NA, kappa3 = NA,
      eff_rank   = NA, sparsity = NA,
      comp_kappa = NA
    ))
  }
  
  list(
    tissue     = tissue,
    samples    = ncol(mat),
    kappa1     = get_cov_kappa(mat),
    kappa2     = get_2norm_kappa(mat),
    kappa3     = get_svd_kappa(mat),
    eff_rank   = get_effective_rank(mat),
    sparsity   = get_matrix_sparsity(mat),
    comp_kappa = compute_composite_kappa(mat)
  )
}
