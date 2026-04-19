#' This script quantifies paired-state transcriptomic entropy using both Shannon
#' and spectral formulations, allowing malignant transformation to be assessed
#' in terms of changes in uncertainty and spectral disorder. In combination with
#' permutation testing, it provides a complementary inferential framework to the
#' complexity engine for evaluating potentially chaotic or anti-chaotic
#' biological processes.
NULL
# -------------------------------------------------------------
# Pairwise Entropy Comparisons
# -------------------------------------------------------------

source(here::here("R/engines/entropy/core_entropy_metrics.R"))
source(here::here("R/engines/entropy/statistical_entropy_helpers.R"))

compare_matrix_pair_entropy <- function(label, chip, mat1, mat2,
                                        filter_probes = NULL,
                                        entropy_fn = compute_shannon_entropy,
                                        n_perm = 1000) {
  if (!is.null(filter_probes)) {
    mat1 <- mat1[rownames(mat1) %in% filter_probes, , drop = FALSE]
    mat2 <- mat2[rownames(mat2) %in% filter_probes, , drop = FALSE]
    
    if (nrow(mat1) < 5 || nrow(mat2) < 5) {
      warning(glue::glue("Too few probes after filtering in {label}"))
      return(NULL)
    }
  }
  
  # Shannon entropy
  shannon_1 <- compute_shannon_entropy(mat1)
  shannon_2 <- compute_shannon_entropy(mat2)
  shannon_delta <- round(shannon_2 - shannon_1, 2)
  shannon_direction <- entropy_direction(shannon_delta)
  
  # Spectral entropy
  spectral_1 <- compute_spectral_entropy(mat1)
  spectral_2 <- compute_spectral_entropy(mat2)
  spectral_delta <- round(spectral_2 - spectral_1, 2)
  spectral_direction <- entropy_direction(spectral_delta)
  
  # Perm test (on Shannon by default)
  perm <- perm_test_entropy(mat1, mat2, entropy_fn = entropy_fn, n_perm = n_perm)
  
  list(
    Comparison = label,
    Chip       = chip,
    Filter     = if (is.null(filter_probes)) "ALL" else paste(length(filter_probes), "probes"),
    Shannon_1  = shannon_1,
    Shannon_2  = shannon_2,
    Shannon_Delta = shannon_delta,
    Shannon_Direction = shannon_direction,
    Spectral_1  = spectral_1,
    Spectral_2  = spectral_2,
    Spectral_Delta = spectral_delta,
    Spectral_Direction = spectral_direction,
    p_perm     = perm$p_value
  )
}
