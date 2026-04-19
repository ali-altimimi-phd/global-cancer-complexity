#' This script quantifies transcriptomic complexity for paired biological states
#' using an SVD-based complexity framework supplemented by permutation testing,
#' bootstrap confidence intervals, per-sample distributional comparisons,
#' effective rank, sparsity, and composite metrics. It thereby operationalizes
#' complexity gain or loss as a multi-criterion structural property of the
#' expression matrix rather than as a single scalar summary.
NULL
# -------------------------------------------------------------
# Pairwise Complexity Comparisons
# -------------------------------------------------------------

source(here::here("R/engines/complexity/core_complexity_metrics.R"))
source(here::here("R/engines/complexity/statistical_complexity_helpers.R"))

compare_matrix_pair_complexity <- function(label, chip, mat1, mat2,
                                           filter_probes = NULL,
                                           complexity_fn = get_svd_kappa,
                                           n_perm = 1000,
                                           n_boot = 1000) {
  # If a filter is provided, subset matrices
  if (!is.null(filter_probes)) {
    mat1 <- mat1[rownames(mat1) %in% filter_probes, , drop = FALSE]
    mat2 <- mat2[rownames(mat2) %in% filter_probes, , drop = FALSE]
    
    # Check for probe dropout
    if (nrow(mat1) < 5 || nrow(mat2) < 5) {
      warning(glue::glue("Too few probes after filtering in {label}"))
      return(NULL)
    }
  }
  
  # Proceed as before
  kappa1 <- complexity_fn(mat1)
  kappa2 <- complexity_fn(mat2)
  delta  <- round(kappa2 - kappa1, 3)
  gain   <- complexity_gain(delta)
  
  perm   <- permutation_test_complexity(mat1, mat2, complexity_fn, n_perm)
  boot_1 <- bootstrap_kappa_ci(mat1, complexity_fn, n_boot)
  boot_2 <- bootstrap_kappa_ci(mat2, complexity_fn, n_boot)
  
  per_sample_1 <- kappa_per_sample(mat1, complexity_fn)
  per_sample_2 <- kappa_per_sample(mat2, complexity_fn)
  
  ks_p    <- suppressWarnings(compare_kappa_distributions(per_sample_1, per_sample_2, "ks"))
  wil_p   <- suppressWarnings(compare_kappa_distributions(per_sample_1, per_sample_2, "wilcox"))
  
  list(
    Comparison = label,
    Chip       = chip,
    Filter     = if (is.null(filter_probes)) "ALL" else paste0(length(filter_probes), " probes"),
    `SVD κ 1`  = kappa1,
    `SVD κ 2`  = kappa2,
    `Δ SVD κ`  = delta,
    Direction  = gain,
    p_perm     = perm$p_value,
    CI_1_Low   = boot_1$ci_lower,
    CI_1_High  = boot_1$ci_upper,
    CI_2_Low   = boot_2$ci_lower,
    CI_2_High  = boot_2$ci_upper,
    p_ks       = ks_p,
    p_wilcox   = wil_p,
    EffRank_1  = get_effective_rank(mat1),
    EffRank_2  = get_effective_rank(mat2),
    Sparsity_1 = get_matrix_sparsity(mat1),
    Sparsity_2 = get_matrix_sparsity(mat2),
    κ_composite_1 = compute_composite_kappa(mat1),
    κ_composite_2 = compute_composite_kappa(mat2),
    perm_dist  = perm$null_distribution,
    boot_1     = boot_1$distribution,
    boot_2     = boot_2$distribution,
    per_sample_1 = per_sample_1,
    per_sample_2 = per_sample_2
  )
}