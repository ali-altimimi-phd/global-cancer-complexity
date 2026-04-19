#' This script defines the core entropy measures used to assess transcriptomic
#' disorder across biological states.
#' 
#' By combining Shannon entropy with spectral entropy derived from covariance 
#' eigenvalues, it treats malignant transformation as a potential shift in both 
#' expression-level uncertainty and global variance organization, and classifies 
#' those shifts as chaotic, anti-chaotic, or neutral.
NULL
# -------------------------------------------------------------
# Core Entropy Metric Functions
# -------------------------------------------------------------
#' Compute Shannon Entropy of a matrix
#'
#' Calculates the Shannon entropy (base 2) of the input numeric matrix,
#' representing the uncertainty or randomness in gene expression data.
#' Returns NA if the input is invalid.
#'
#' @param mat Numeric matrix of gene expression values
#' @return Numeric Shannon entropy value
#' @export
compute_shannon_entropy <- function(mat) {
  if (!is.matrix(mat) || nrow(mat) == 0 || ncol(mat) == 0) return(NA_real_)
  mat <- na.omit(mat)
  round(DescTools::Entropy(mat, base = 2), 2)
}

#' Compute Spectral Entropy of a matrix
#'
#' Calculates spectral entropy based on eigenvalues of the covariance matrix,
#' capturing signal complexity in terms of spectral content and variance distribution.
#' Returns NA if input is invalid or covariance fails.
#'
#' @param mat Numeric matrix of gene expression values
#' @return Numeric spectral entropy value
#' @export
compute_spectral_entropy <- function(mat) { ... }

#' Classify Entropy Direction from difference values
#'
#' Translates entropy change (delta) into descriptive categories:
#' strongly anti-chaotic, mildly anti-chaotic, neutral, mildly chaotic, strongly chaotic.
#' Useful for biological interpretation of tumor state shifts.
#'
#' @param delta Numeric difference in entropy between two conditions
#' @return Character entropy category label
#' @export
entropy_direction <- function(delta) { ... }

# -------------------------------------------------------------
# Core Entropy Metric Functions
# -------------------------------------------------------------

library(dplyr)

compute_spectral_entropy <- function(mat) {
  if (!is.matrix(mat) || nrow(mat) == 0 || ncol(mat) == 0) return(NA_real_)
  mat <- na.omit(mat)
  
  cov_mat <- tryCatch(cov(mat), error = function(e) return(NULL))
  eigs <- tryCatch(eigen(cov_mat)$values, error = function(e) return(NULL))
  if (is.null(eigs)) return(NA_real_)
  
  eigs <- eigs / sum(eigs + .Machine$double.eps)
  entropy <- -sum(eigs * log2(eigs + .Machine$double.eps))
  round(entropy, 3)
}

entropy_direction <- function(delta) {
  case_when(
    is.na(delta) ~ NA_character_,
    delta < -0.5 ~ "strongly anti-chaotic",
    delta < -0.1 ~ "mildly anti-chaotic",
    delta <  0.1 ~ "neutral",
    delta <  0.5 ~ "mildly chaotic",
    TRUE         ~ "strongly chaotic"
  )
}
