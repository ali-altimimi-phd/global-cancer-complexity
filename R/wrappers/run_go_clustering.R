# ------------------------------------------------------------------------------
# File: run_go_clustering.R
# Purpose: Orchestrate GO clustering analysis across all comparisons and manage
#   execution of comparison-level clustering workflows
# Role: Analysis orchestration utility
# Pipeline: Reporting
# Project: Cancer Complexity Analysis
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------

#' Run GO Clustering for All Comparisons
#'
#' Iterates over all comparisons in `summaries_combined_df` and runs GO
#' clustering analysis using semantic similarity–based grouping.
#'
#' @param summaries_combined_df Data frame containing comparison-level summaries.
#'   The `comparison` column is used to determine which comparisons to process.
#' @param complexity_df Data frame. Cleaned complexity results (GO, KEGG, etc.).
#' @param entropy_df Data frame. Cleaned entropy results (reserved for future use).
#' @param output_dir_base Character. Base directory where HTML summaries are written.
#' @param go_mode Character. Optional GO ontology filter. One of
#'   `"GO_BP"`, `"GO_MF"`, or `"GO_UNSPECIFIED"`. Passed through to
#'   downstream clustering functions.
run_go_clustering_main <- function(summaries_combined_df,
                                   complexity_df,
                                   entropy_df,
                                   output_dir_base = "quarto/resources/tables/clusters",
                                   go_mode = NULL) {
  
  source(here::here("R/analyze/analyze_comparison_terms.R"))
  source(here::here("R/helpers/sanitize_comparison.R"))
  
  all_comparisons <- unique(summaries_combined_df$comparison)

  for (cmp in all_comparisons) {
    message("→ Generating GO analysis for: ", cmp)
    
    try({
      analyze_comparison_terms(
        cmp = cmp,
        mode = "GO",
        go_mode = go_mode,
        complexity_df = complexity_df,
        entropy_df = entropy_df,
        output_dir_base = output_dir_base
      )
    }, silent = TRUE)
  }
}
