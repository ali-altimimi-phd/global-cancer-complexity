# -------------------------------------------------------------
# Matrix Comparison Across Map
# -------------------------------------------------------------

run_pairwise_entropy <- function(comparison_df,
                                 entropy_fn = compute_shannon_entropy,
                                 n_perm = 1000) {
  stopifnot(all(c("comparison_label", "chip", "matrix_1", "matrix_2") %in% colnames(comparison_df)))
  
  purrr::pmap_dfr(comparison_df, function(comparison_label, chip, matrix_1, matrix_2, filter_probes = NULL) {
    result <- compare_matrix_pair_entropy(
      label = comparison_label,
      chip = chip,
      mat1 = matrix_1,
      mat2 = matrix_2,
      filter_probes = filter_probes,
      entropy_fn = entropy_fn,
      n_perm = n_perm
    )
    
    if (!is.null(result)) tibble::as_tibble(result) else NULL
  })
}


