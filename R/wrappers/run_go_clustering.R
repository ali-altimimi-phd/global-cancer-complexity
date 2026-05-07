# ------------------------------------------------------------------------------
# File: run_go_clustering.R
# Purpose: Orchestrate ontology-aware GO semantic clustering for reporting.
# Role: Reporting pipeline wrapper
# Pipeline: Reporting
# Project: Cancer Complexity Analysis
# Author: Ali M. Al-Timimi; revised with ChatGPT assistance
# Created: 2026
# ------------------------------------------------------------------------------

#' Run GO Semantic Clustering for Reporting
#'
#' This wrapper replaces the old per-comparison dispatcher. It runs the unified
#' GO semantic summarization layer once across all requested comparisons, writes
#' aggregate outputs, and also writes per-comparison RDS/CSV files for Quarto.
#'
#' @param summaries_combined_df Data frame containing a `comparison` column.
#' @param complexity_df Data frame. Cleaned complexity results.
#' @param entropy_df Data frame. Cleaned entropy results.
#' @param output_dir_base Character. Base directory where structured outputs are written.
#' @param go_mode Character. Optional GO ontology filter: GO_ALL, GO_BP, GO_MF, GO_CC.
#' @param alpha Numeric. Term-level p_perm cutoff.
#' @param similarity_cutoff Numeric. GO semantic similarity clustering cutoff.
#' @param write_html Logical. Optional convenience HTML table output.
#' @param logger Optional pipeline logger with a $log method.
#'
#' @return Invisibly returns a list with significant_go_terms, clustered_terms,
#'   and cluster_summary.
run_go_clustering_main <- function(summaries_combined_df,
                                   complexity_df,
                                   entropy_df,
                                   output_dir_base = here::here("quarto", "resources", "tables", "clusters"),
                                   go_mode = NULL,
                                   alpha = 0.05,
                                   similarity_cutoff = 0.70,
                                   write_html = FALSE,
                                   logger = NULL) {
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package 'here' is required by the reporting pipeline wrapper.", call. = FALSE)
  }

  source(here::here("R", "helpers", "sanitize_comparison.R"))
  source(here::here("R", "analyze_go_clusters", "analyze_go_semantic_clusters.R"))

  log_msg <- function(...) {
    msg <- paste0(...)
    if (!is.null(logger) && is.function(logger$log)) {
      logger$log(msg, section = "GO_CLUSTER")
    } else {
      message(msg)
    }
  }

  if (!"comparison" %in% names(summaries_combined_df)) {
    stop("summaries_combined_df must contain a 'comparison' column.", call. = FALSE)
  }

  all_comparisons <- unique(stats::na.omit(as.character(summaries_combined_df$comparison)))

  if (length(all_comparisons) == 0) {
    stop("No comparisons found in summaries_combined_df$comparison.", call. = FALSE)
  }

  log_msg("Running GO semantic clustering for ", length(all_comparisons), " comparisons.")
  log_msg("GO mode: ", ifelse(is.null(go_mode), "GO_ALL", go_mode))
  log_msg("Term-level alpha: ", alpha)
  log_msg("Semantic similarity cutoff: ", similarity_cutoff)

  res <- analyze_go_semantic_clusters(
    complexity_df = complexity_df,
    entropy_df = entropy_df,
    comparisons = all_comparisons,
    output_dir = output_dir_base,
    alpha = alpha,
    go_mode = go_mode,
    similarity_cutoff = similarity_cutoff,
    write_html = write_html,
    logger = logger
  )

  log_msg("GO semantic clustering complete.")
  log_msg("Significant GO rows: ", nrow(res$significant_go_terms))
  log_msg("Clustered GO rows: ", nrow(res$clustered_terms))
  log_msg("Cluster summary rows: ", nrow(res$cluster_summary))

  invisible(res)
}
