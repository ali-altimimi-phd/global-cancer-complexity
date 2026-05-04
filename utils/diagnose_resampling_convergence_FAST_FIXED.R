# ------------------------------------------------------------------------------
# File: diagnose_resampling_convergence.R
# Purpose: Estimate convergence of resampled complexity/entropy metrics across
#   chip x comparison FULL limma/variance-filtered probe sets.
# Role: Standalone diagnostic utility
# Pipeline: Analysis diagnostics
# Project: Cancer Complexity Analysis
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------

#' Diagnose Resampling Convergence for Pairwise Cancer Comparisons
#'
#' @description
#' This standalone utility reads the filtered probe objects produced by the
#' comparison-aware filtering stage and evaluates how eigenstructure-based
#' complexity/entropy summaries stabilize as the number of resampling iterations
#' increases.
#'
#' The diagnostic is intentionally read-only with respect to existing pipeline
#' outputs. It writes only diagnostic CSV files to
#' `output/<study_name>/diagnostics/resampling/`.
#'
#' @details
#' Diagnostic unit:
#'   chip x comparison x filtered FULL probe set
#'
#' Metrics:
#'   - participation ratio: complexity / effective dimensionality proxy
#'   - eigenvalue entropy: entropy / spectral dispersion proxy
#'
#' Resampling definition:
#'   - probes are resampled without replacement from the comparison-specific
#'     filtered probe set
#'   - each resampling level is run independently
#'
#' @note
#' This script is designed to be tolerant of moderate object-structure changes
#' in `filtered_probes_*.rds`, but the helper `extract_diagnostic_matrix()` may
#' need light adaptation if your saved filtering object uses different names.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

# ---- Configuration + logging ----
source(here::here("R/config/global_cancer/analysis_config.R"), local = FALSE)
source(here::here("R/helpers/pipeline_logger.R"))

resampling_logfile <- file.path(
  analysis_logs_dir,
  "resampling_convergence_diagnostic_log.txt"
)

logger <- start_log(resampling_logfile)
logger$log("Starting resampling convergence diagnostic.", section = "RESAMPLING")

# ---- Diagnostic parameters ----
# Keep this grid small enough to run quickly, but broad enough to evaluate whether
# the historical/default value of 1000 is necessary.
resample_grid <- c(50, 100, 250)

# Probe fraction used inside each resample. A fraction of 0.80 means each replicate
# recomputes metrics on 80% of the filtered FULL probe set for that comparison.
probe_fraction <- 0.80

# Optional lower bound. Comparisons below this threshold are skipped because the
# covariance/eigenstructure estimates become difficult to interpret.
min_filtered_probes_for_diagnostic <- 10

# Base seed for reproducibility. Each chip x comparison x resampling-level receives
# a deterministic derived seed.
base_seed <- 20260501L

# Output directory
resampling_diagnostic_dir <- here::here(
  "output", study_name, "diagnostics", "resampling"
)
dir.create(resampling_diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Metric helpers ----
participation_ratio <- function(eigvals) {
  eigvals <- eigvals[is.finite(eigvals) & eigvals > 0]
  if (length(eigvals) == 0 || sum(eigvals) <= 0) return(NA_real_)
  (sum(eigvals)^2) / sum(eigvals^2)
}

eig_entropy <- function(eigvals) {
  eigvals <- eigvals[is.finite(eigvals) & eigvals > 0]
  if (length(eigvals) == 0 || sum(eigvals) <= 0) return(NA_real_)
  p <- eigvals / sum(eigvals)
  -sum(p * log(p))
}

compute_spectral_metrics <- function(expr_mat) {
  # Expected orientation: probes x samples.
  # We compute covariance among samples across selected probes.
  # IMPORTANT: for a probes x samples matrix, cov(expr_mat) is samples x samples.
  # Do NOT use cov(t(expr_mat)); that creates a probes x probes matrix and can hang.
  if (!is.matrix(expr_mat)) expr_mat <- as.matrix(expr_mat)

  if (nrow(expr_mat) < 2 || ncol(expr_mat) < 2) {
    return(tibble(
      participation_ratio = NA_real_,
      eig_entropy = NA_real_
    ))
  }

  cov_mat <- stats::cov(expr_mat, use = "pairwise.complete.obs")
  eigvals <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
  eigvals <- pmax(eigvals, 0)

  tibble(
    participation_ratio = participation_ratio(eigvals),
    eig_entropy = eig_entropy(eigvals)
  )
}

# ---- Object-structure helpers ----
`%||%` <- function(x, y) if (!is.null(x)) x else y

is_expression_matrix <- function(x) {
  is.matrix(x) && is.numeric(x) && nrow(x) >= 2 && ncol(x) >= 2
}

find_first_matrix <- function(x, max_depth = 4) {
  if (max_depth < 0) return(NULL)
  if (is_expression_matrix(x)) return(x)

  if (is.data.frame(x)) {
    numeric_cols <- vapply(x, is.numeric, logical(1))
    if (sum(numeric_cols) >= 2 && nrow(x) >= 2) {
      candidate <- as.matrix(x[, numeric_cols, drop = FALSE])
      if (is_expression_matrix(candidate)) return(candidate)
    }
  }

  if (is.list(x)) {
    # Prefer likely expression-matrix field names first.
    preferred_names <- c(
      "expr_mat", "expr_matrix", "expression_matrix", "filtered_matrix",
      "matrix", "data", "x", "filtered_expr", "filtered_expression"
    )

    nm <- names(x)
    if (!is.null(nm)) {
      for (candidate_name in intersect(preferred_names, nm)) {
        found <- find_first_matrix(x[[candidate_name]], max_depth = max_depth - 1)
        if (!is.null(found)) return(found)
      }
    }

    for (item in x) {
      found <- find_first_matrix(item, max_depth = max_depth - 1)
      if (!is.null(found)) return(found)
    }
  }

  NULL
}

extract_selected_probes <- function(x) {
  if (!is.list(x)) return(NULL)

  preferred_names <- c(
    "selected_probes", "filtered_probes", "probe_ids", "probes",
    "probe_id", "ids", "features", "feature_ids"
  )

  nm <- names(x)
  if (is.null(nm)) return(NULL)

  for (candidate_name in intersect(preferred_names, nm)) {
    candidate <- x[[candidate_name]]
    if (is.character(candidate) || is.factor(candidate)) {
      candidate <- as.character(candidate)
      if (length(candidate) > 0) return(candidate)
    }
    if (is.data.frame(candidate)) {
      for (col in names(candidate)) {
        if (is.character(candidate[[col]]) || is.factor(candidate[[col]])) {
          vals <- as.character(candidate[[col]])
          if (length(vals) > 0) return(vals)
        }
      }
    }
  }

  NULL
}

extract_sample_counts <- function(x, expr_mat = NULL) {
  out <- list(n_samples_normal = NA_integer_, n_samples_tumor = NA_integer_)

  if (is.list(x)) {
    nm <- names(x)
    if (!is.null(nm)) {
      normal_names <- c("normal_samples", "normal_sample_ids", "normal_cols", "normal")
      tumor_names  <- c("tumor_samples", "tumor_sample_ids", "tumor_cols", "tumor")

      for (candidate_name in intersect(normal_names, nm)) {
        candidate <- x[[candidate_name]]
        if (is.vector(candidate) || is.list(candidate)) {
          out$n_samples_normal <- length(candidate)
          break
        }
      }
      for (candidate_name in intersect(tumor_names, nm)) {
        candidate <- x[[candidate_name]]
        if (is.vector(candidate) || is.list(candidate)) {
          out$n_samples_tumor <- length(candidate)
          break
        }
      }
    }
  }

  if (is.na(out$n_samples_normal) && !is.null(expr_mat)) {
    out$n_samples_normal <- NA_integer_
  }
  if (is.na(out$n_samples_tumor) && !is.null(expr_mat)) {
    out$n_samples_tumor <- NA_integer_
  }

  out
}

extract_diagnostic_matrix <- function(comparison_obj) {
  expr_mat <- find_first_matrix(comparison_obj)
  if (is.null(expr_mat)) return(NULL)

  selected_probes <- extract_selected_probes(comparison_obj)

  # If selected probes exist and matrix rownames are available, enforce that the
  # diagnostic uses only the filtered FULL probe set.
  if (!is.null(selected_probes) && !is.null(rownames(expr_mat))) {
    keep <- intersect(rownames(expr_mat), selected_probes)
    if (length(keep) >= 2) {
      expr_mat <- expr_mat[keep, , drop = FALSE]
    }
  }

  storage.mode(expr_mat) <- "double"
  expr_mat
}

# ---- Input loading ----
load_filtered_probe_object <- function(chip_id) {
  rds_path <- file.path(filtered_probes_dir, sprintf("filtered_probes_%s.rds", chip_id))

  if (!file.exists(rds_path)) {
    stop(sprintf(
      "Missing filtered-probe object for %s: %s. Run filtering first.",
      chip_id, rds_path
    ), call. = FALSE)
  }

  readRDS(rds_path)
}

as_comparison_list <- function(filtered_obj) {
  # Most likely case: top-level list where each element is a comparison.
  if (is.list(filtered_obj) && !is.null(names(filtered_obj))) {
    # If the object has an obvious nested comparison list, prefer it.
    for (candidate_name in c("comparisons", "results", "filtered_results", "comparison_results")) {
      if (candidate_name %in% names(filtered_obj) && is.list(filtered_obj[[candidate_name]])) {
        return(filtered_obj[[candidate_name]])
      }
    }
    return(filtered_obj)
  }

  stop("Filtered object is not a named list; inspect filtered_probes_*.rds structure.", call. = FALSE)
}

# ---- Resampling core ----
run_one_resampling_level <- function(expr_mat,
                                     n_resamples,
                                     probe_fraction,
                                     seed) {
  set.seed(seed)

  n_probes_full <- nrow(expr_mat)
  n_probes_sampled <- max(2L, floor(n_probes_full * probe_fraction))
  n_probes_sampled <- min(n_probes_sampled, n_probes_full)

  start_time <- Sys.time()

  metric_tbl <- map_dfr(seq_len(n_resamples), function(i) {
    sampled_idx <- sample.int(
      n = n_probes_full,
      size = n_probes_sampled,
      replace = FALSE
    )

    compute_spectral_metrics(expr_mat[sampled_idx, , drop = FALSE]) |>
      mutate(iteration = i, .before = 1)
  })

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  metric_tbl |>
    summarise(
      n_resamples = n_resamples,
      n_probes_sampled = n_probes_sampled,
      mean_pr = mean(participation_ratio, na.rm = TRUE),
      sd_pr = sd(participation_ratio, na.rm = TRUE),
      cv_pr = sd_pr / mean_pr,
      mean_entropy = mean(eig_entropy, na.rm = TRUE),
      sd_entropy = sd(eig_entropy, na.rm = TRUE),
      cv_entropy = sd_entropy / mean_entropy,
      runtime_seconds = elapsed,
      .groups = "drop"
    )
}

run_comparison_diagnostic <- function(chip_id,
                                      comparison_name,
                                      comparison_obj,
                                      comparison_index) {
  expr_mat <- extract_diagnostic_matrix(comparison_obj)

  if (is.null(expr_mat)) {
    logger$log(
      sprintf("Skipping %s / %s: no expression matrix found.", chip_id, comparison_name),
      section = "WARNING"
    )
    return(NULL)
  }

  n_probes_full <- nrow(expr_mat)
  n_samples_total <- ncol(expr_mat)
  sample_counts <- extract_sample_counts(comparison_obj, expr_mat)

  if (n_probes_full < min_filtered_probes_for_diagnostic) {
    logger$log(
      sprintf(
        "Skipping %s / %s: only %d filtered probes.",
        chip_id, comparison_name, n_probes_full
      ),
      section = "WARNING"
    )
    return(NULL)
  }

  logger$log(
    sprintf(
      "Running %s / %s: %d probes, %d samples.",
      chip_id, comparison_name, n_probes_full, n_samples_total
    ),
    section = "RESAMPLING"
  )

  map_dfr(seq_along(resample_grid), function(grid_index) {
    n_resamples <- resample_grid[[grid_index]]

    logger$log(
      sprintf(
        "  Resampling level: %s / %s | n_resamples = %d",
        chip_id, comparison_name, n_resamples
      ),
      section = "RESAMPLING"
    )

    derived_seed <- base_seed +
      (comparison_index * 1000L) +
      (grid_index * 10L) +
      match(chip_id, chips)

    run_one_resampling_level(
      expr_mat = expr_mat,
      n_resamples = n_resamples,
      probe_fraction = probe_fraction,
      seed = derived_seed
    ) |>
      mutate(
        chip = chip_id,
        comparison = comparison_name,
        n_samples_total = n_samples_total,
        n_samples_normal = sample_counts$n_samples_normal,
        n_samples_tumor = sample_counts$n_samples_tumor,
        n_probes_full = n_probes_full,
        probe_fraction = probe_fraction,
        seed = derived_seed,
        .before = 1
      )
  })
}

# ---- Run diagnostics ----
all_results <- list()

for (chip_id in chips) {
  logger$log(sprintf("Loading filtered object for chip %s.", chip_id), section = "RESAMPLING")

  filtered_obj <- load_filtered_probe_object(chip_id)
  comparison_list <- as_comparison_list(filtered_obj)

  if (length(comparison_list) == 0) {
    logger$log(sprintf("No comparisons found for chip %s.", chip_id), section = "WARNING")
    next
  }

  comparison_names <- names(comparison_list)
  if (is.null(comparison_names)) {
    comparison_names <- paste0("comparison_", seq_along(comparison_list))
  }

  chip_results <- map2_dfr(
    comparison_list,
    seq_along(comparison_list),
    function(comparison_obj, comparison_index) {
      comparison_name <- comparison_names[[comparison_index]]
      run_comparison_diagnostic(
        chip_id = chip_id,
        comparison_name = comparison_name,
        comparison_obj = comparison_obj,
        comparison_index = comparison_index
      )
    }
  )

  all_results[[chip_id]] <- chip_results
}

resampling_detail_tbl <- bind_rows(all_results)

if (nrow(resampling_detail_tbl) == 0) {
  stop("No diagnostic results were produced. Inspect filtered object structure.", call. = FALSE)
}

# ---- Summaries ----
resampling_global_summary_tbl <- resampling_detail_tbl |>
  group_by(n_resamples) |>
  summarise(
    n_runs = n(),
    median_cv_pr = median(cv_pr, na.rm = TRUE),
    q75_cv_pr = quantile(cv_pr, probs = 0.75, na.rm = TRUE),
    max_cv_pr = max(cv_pr, na.rm = TRUE),
    median_cv_entropy = median(cv_entropy, na.rm = TRUE),
    q75_cv_entropy = quantile(cv_entropy, probs = 0.75, na.rm = TRUE),
    max_cv_entropy = max(cv_entropy, na.rm = TRUE),
    total_runtime_seconds = sum(runtime_seconds, na.rm = TRUE),
    .groups = "drop"
  )

# Loose diagnostic thresholds; these do not change the pipeline. They simply help
# choose a defensible default for the config later.
cv_threshold_pr <- 0.05
cv_threshold_entropy <- 0.05

resampling_stability_tbl <- resampling_detail_tbl |>
  mutate(
    stable_pr = cv_pr <= cv_threshold_pr,
    stable_entropy = cv_entropy <= cv_threshold_entropy,
    stable_both = stable_pr & stable_entropy
  ) |>
  group_by(n_resamples) |>
  summarise(
    n_runs = n(),
    n_stable_pr = sum(stable_pr, na.rm = TRUE),
    n_stable_entropy = sum(stable_entropy, na.rm = TRUE),
    n_stable_both = sum(stable_both, na.rm = TRUE),
    pct_stable_both = 100 * n_stable_both / n_runs,
    .groups = "drop"
  )

# ---- Write outputs ----
detail_path <- file.path(
  resampling_diagnostic_dir,
  "resampling_convergence_by_comparison.csv"
)
summary_path <- file.path(
  resampling_diagnostic_dir,
  "resampling_convergence_global_summary.csv"
)
stability_path <- file.path(
  resampling_diagnostic_dir,
  "resampling_convergence_stability_summary.csv"
)

write_csv(resampling_detail_tbl, detail_path)
write_csv(resampling_global_summary_tbl, summary_path)
write_csv(resampling_stability_tbl, stability_path)

logger$log(sprintf("Wrote comparison-level diagnostics: %s", detail_path), section = "RESAMPLING")
logger$log(sprintf("Wrote global summary: %s", summary_path), section = "RESAMPLING")
logger$log(sprintf("Wrote stability summary: %s", stability_path), section = "RESAMPLING")
logger$log("Completed resampling convergence diagnostic.", section = "RESAMPLING")

# Print compact summary for interactive/RStudio use.
print(resampling_global_summary_tbl)
print(resampling_stability_tbl)
