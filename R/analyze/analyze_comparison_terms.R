# ------------------------------------------------------------------------------
# File: analyze_comparison_terms.R
# Purpose: Dispatch comparison-level gene set analysis across GO, KEGG, or MSigDB modes
# Role: Analysis dispatch utility
# Pipeline: Reporting
# Project: Cancer Complexity Analysis
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------

#' Analyze gene set clustering for a comparison
#'
#' This function runs clustering and summarization of GO, KEGG, or MSigDB gene set results
#' for a specific cancer comparison and saves the resulting summary tables to disk.
#'
#' Currently, only GO clustering is implemented. KEGG and MSIG modes are placeholders.
#'
#' @param cmp Character. The comparison label (e.g., `"BR/BRAD"`).
#' @param mode Character. Gene set category to analyze. Must be one of `"GO"`, `"KEGG"`, or `"MSIG"`.
#' @param complexity_df Data frame. Complexity results for all comparisons and sets.
#' @param entropy_df Data frame. Entropy results for all comparisons and sets.
#' @param output_dir_base Character. Path to the base directory for saving output tables.
#' @param go_mode Character. Optional GO ontology filter: one of `"GO_BP"`, `"GO_MF"`, or `"GO_UNSPECIFIED"`.
#'
#' @return Invisibly returns `NULL`. Output tables are written to disk as side effects.
#' 
#' @seealso [analyze_go_complexity_clusters()] for GO complexity clustering logic.
analyze_comparison_terms <- function(cmp,
                                     mode,
                                     complexity_df,
                                     entropy_df,
                                     output_dir_base = "quarto/resources/tables/clusters",
                                     go_mode = NULL) {
  # Create subdir per mode
  output_dir <- file.path(output_dir_base, tolower(mode))
  fs::dir_create(output_dir)
  
  source(here::here("R/analyze/analyze_go_complexity_clusters.R"))  # renamed file
  # source(here::here("R/analyze/analyze_go_entropy_clusters.R"))
  
  # Dispatch by mode
  if (mode == "GO") {
    analyze_go_complexity_clusters(
      cmp = cmp,
      complexity_df = complexity_df,
      output_dir = output_dir,
      go_mode = go_mode
    )
    
    # analyze_go_entropy_clusters(
    #   cmp = cmp,
    #   entropy_df = entropy_df,
    #   output_dir = output_dir,
    #   go_mode = go_mode
    # )
  } else if (mode == "KEGG") {
    message("KEGG clustering not implemented yet.")
  } else if (mode == "MSIG") {
    message("MSIG analysis not implemented yet.")
  } else {
    stop("Unsupported mode: ", mode)
  }
  
  invisible(NULL)
}
