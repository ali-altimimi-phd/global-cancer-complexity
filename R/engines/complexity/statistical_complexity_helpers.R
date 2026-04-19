#' This script provides inferential support for matrix-level complexity analysis
#' through bootstrap confidence intervals, permutation testing, and distributional
#' comparison utilities. These procedures allow observed differences in transcriptomic
#' complexity between biological states to be evaluated for stability and statistical
#' plausibility rather than treated as raw descriptive contrasts alone.
NULL
# -------------------------------------------------------------
# Statistical Helper Functions
# -------------------------------------------------------------

#' Bootstrap confidence interval for a complexity metric
#'
#' @param mat A numeric matrix.
#' @param complexity_fn Function to compute a complexity score (default: get_svd_kappa).
#' @param n_boot Number of bootstrap samples (default: 1000).
#' @param level Confidence level (default: 0.95).
#' @param sample_dim Whether to resample columns ("col") or rows ("row").
#' @return A list with mean, lower and upper confidence bounds, and the bootstrap distribution.
#' @export
bootstrap_kappa_ci <- function(mat, complexity_fn = get_svd_kappa,
                               n_boot = 1000, level = 0.95, sample_dim = "col") {
  stats <- replicate(n_boot, {
    if (sample_dim == "col") {
      mat_boot <- mat[, sample(ncol(mat), replace = TRUE)]
    } else {
      mat_boot <- mat[sample(nrow(mat), replace = TRUE), ]
    }
    complexity_fn(mat_boot)
  })
  
  ci <- quantile(stats, probs = c((1 - level) / 2, 1 - (1 - level) / 2))
  
  list(
    mean = mean(stats),
    ci_lower = ci[[1]],
    ci_upper = ci[[2]],
    distribution = stats
  )
}

#' Compute complexity score per sample
#'
#' @param mat A numeric matrix.
#' @param complexity_fn Function to apply to each column (default: get_svd_kappa).
#' @return Numeric vector of complexity scores for each sample.
#' @export
kappa_per_sample <- function(mat, complexity_fn = get_svd_kappa) {
  if (is.null(mat) || ncol(mat) < 2) return(numeric(0))
  apply(mat, 2, function(col) complexity_fn(matrix(col, nrow = length(col))))
}

#' Compare distributions of complexity scores
#'
#' @param norm_vals Vector of normal sample scores.
#' @param cancer_vals Vector of cancer sample scores.
#' @param method "ks" for Kolmogorov-Smirnov, "wilcox" for Wilcoxon test.
#' @return P-value from the selected statistical test.
#' @export
compare_kappa_distributions <- function(norm_vals, cancer_vals, method = "ks") {
  if (method == "ks") {
    test <- ks.test(norm_vals, cancer_vals)
  } else {
    test <- wilcox.test(norm_vals, cancer_vals)
  }
  test$p.value
}

#' Permutation test for complexity score difference
#'
#' @param mat_normal Matrix of normal samples.
#' @param mat_cancer Matrix of cancer samples.
#' @param complexity_fn Complexity function to apply (default: get_svd_kappa).
#' @param n_perm Number of permutations (default: 1000).
#' @return A list with observed difference, p-value, and null distribution.
#' @export
permutation_test_complexity <- function(mat_normal, mat_cancer,
                                        complexity_fn = get_svd_kappa,
                                        n_perm = 1000) {
  all_mat <- cbind(mat_normal, mat_cancer)
  labels <- c(rep("normal", ncol(mat_normal)), rep("cancer", ncol(mat_cancer)))
  
  observed_diff <- complexity_fn(mat_cancer) - complexity_fn(mat_normal)
  
  null_dist <- replicate(n_perm, {
    perm_labels <- sample(labels)
    mat_n <- all_mat[, perm_labels == "normal"]
    mat_c <- all_mat[, perm_labels == "cancer"]
    complexity_fn(mat_c) - complexity_fn(mat_n)
  })
  
  p_value <- mean(abs(null_dist) >= abs(observed_diff))
  
  list(
    observed = observed_diff,
    p_value = p_value,
    null_distribution = null_dist
  )
}
