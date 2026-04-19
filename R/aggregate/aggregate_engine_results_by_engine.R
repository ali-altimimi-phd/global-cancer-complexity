# ------------------------------------------------------------------------------
# File: aggregate_engine_results_by_engine.R
# Purpose: Aggregate comparison-level complexity or entropy results across chips,
#   biological comparisons, and gene-set scopes into a unified tabular structure
# Role: Analysis aggregation utility
# Pipeline: Reporting
# Project: Cancer Complexity Analysis
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------

#' Aggregate Engine Results by Engine
#'
#' Aggregates complexity or entropy `.rds` outputs from the analysis pipeline
#' across chip platforms, biological comparisons, and gene-set scopes into a
#' unified tabular representation.
#'
#' This function flattens nested result objects into a common analytical frame,
#' enabling cross-cancer and cross-process interpretation of transcriptomic
#' complexity and entropy patterns.
#'
#' Large list-like columns may be dropped by default to keep the returned table
#' compact and easier to inspect downstream. Gene set names are not attached
#' here; use `attach_gene_set_names()` afterward if needed.
#'
#' @param engine Character scalar. Either `"complexity"` or `"entropy"`.
#' @param input_dir Character scalar. Directory containing RDS result files.
#'   Defaults to `output/global_cancer/RData`.
#' @param output_dir Character scalar. Directory to write outputs if needed.
#'   Included for pipeline consistency; not used directly here.
#' @param drop_large_columns Logical. If `TRUE` (default), remove bulky columns
#'   such as permutation or bootstrap distributions.
#'
#' @return A tibble combining all detected results with metadata columns such as
#'   `chip`, `engine`, `mode`, and `gene_set`.
aggregate_engine_results_by_engine <- function(engine = c("complexity", "entropy"),
                                               input_dir = here::here("output", "global_cancer", "RData"),
                                               output_dir = here::here(data_dir, "aggregated"),
                                               drop_large_columns = TRUE) {
  engine <- match.arg(engine)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  pattern <- glue::glue("{engine}_results_.*\\.rds$")
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  
  drop_cols <- c("perm_dist", "boot_1", "boot_2", "per_sample_1", "per_sample_2")
  
  safe_extract <- function(x) {
    if (drop_large_columns) {
      x[setdiff(names(x), drop_cols)]
    } else {
      x
    }
  }
  
  all_results <- purrr::map_dfr(files, function(file_path) {
    file_name <- fs::path_file(file_path)
    chip <- stringr::str_extract(file_name, "hu[0-9a-z]+")
    mode <- dplyr::case_when(
      stringr::str_detect(file_name, "all")    ~ "ALL",
      stringr::str_detect(file_name, "go")     ~ "GO",
      stringr::str_detect(file_name, "kegg")   ~ "KEGG",
      stringr::str_detect(file_name, "msigdb") ~ "MSIG",
      TRUE ~ "unknown"
    )
    
    results <- readRDS(file_path)
    
    purrr::imap_dfr(results, function(res_by_set, comparison) {
      purrr::imap_dfr(res_by_set, function(entry, gene_set_id) {
        if (is.null(entry)) return(NULL)
        
        tibble::as_tibble(safe_extract(entry)) |>
          dplyr::mutate(
            chip     = chip,
            engine   = engine,
            mode     = mode,
            gene_set = gene_set_id
          )
      })
    })
  })
  
  return(all_results)
}
