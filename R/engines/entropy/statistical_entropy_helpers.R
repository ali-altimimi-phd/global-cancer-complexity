#' This script provides nonparametric inference for entropy-based comparisons 
#' by generating permutation null distributions for observed state-to-state 
#' entropy differences. 
#' 
#' It thereby allows apparent shifts in transcriptomic disorder to be evaluated 
#' against randomized label assignments rather than interpreted as purely
#' descriptive contrasts.
NULL
# -------------------------------------------------------------
# Statistical Helper Functions for Entropy Analysis
# -------------------------------------------------------------
#' Permutation test for entropy differences
#'
#' Performs a permutation test comparing entropy between two groups.
#' Generates a null distribution by random label shuffling,
#' estimating p-values for observed entropy shifts.
#'
#' @param mat_1 Numeric matrix for group 1
#' @param mat_2 Numeric matrix for group 2
#' @param entropy_fn Function to compute entropy (default: Shannon entropy)
#' @param n_perm Number of permutations for null distribution (default: 1000)
#' @return List with delta (entropy difference), p_value, and null_distribution vector
#' @export
perm_test_entropy <- function(mat_1, mat_2,
                              entropy_fn = compute_shannon_entropy,
                              n_perm = 1000) {
  if (!is.matrix(mat_1) || !is.matrix(mat_2)) {
    return(list(delta = NA, p_value = NA, null_distribution = NA))
  }
  
  all_mat <- cbind(mat_1, mat_2)
  labels <- c(rep("A", ncol(mat_1)), rep("B", ncol(mat_2)))
  
  observed_diff <- entropy_fn(mat_2) - entropy_fn(mat_1)
  
  null_dist <- replicate(n_perm, {
    perm_labels <- sample(labels)
    mat_a <- all_mat[, perm_labels == "A", drop = FALSE]
    mat_b <- all_mat[, perm_labels == "B", drop = FALSE]
    entropy_fn(mat_b) - entropy_fn(mat_a)
  })
  
  p_value <- mean(abs(null_dist) >= abs(observed_diff))
  
  list(
    delta = observed_diff,
    p_value = round(p_value, 4),
    null_distribution = null_dist
  )
}
