# ------------------------------------------------------------------------------
# File: run_grouped_probe_filtering.R
# Purpose: Apply comparison-aware probe selection with provenance metadata
# Role: Analysis wrapper
# Pipeline: Analysis
# Project: Cancer Complexity Analysis
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------

#' Comparison-aware probe selection wrapper
#'
#' Performs comparison-specific reduction of expression space by identifying probes
#' most relevant to each tissue-specific normal-to-tumor transition.
#'
#' The selected probes define the molecular feature space over which downstream
#' complexity and entropy measures are evaluated. This wrapper supports both
#' limma-based differential expression filtering and variance-based filtering.
#'
#' In addition to returning per-comparison filtered matrices and probe sets, the
#' function records cross-comparison probe reuse summaries and global provenance
#' metadata describing the filtering method and parameter values used to generate
#' the object.
#'
#' @param matrix_list A named list of expression matrices from `build_matrix_lists_by_tissue()`.
#' @param comparison_map A named list of group -> label -> tissue comparisons.
#' @param chip_id A string identifying the chip, e.g. `"hu35ksuba"`.
#' @param method Filtering method; either `"limma"` or `"variance"`.
#' @param logfc_cutoff Log2 fold-change threshold used when `method = "limma"`.
#' @param pval_cutoff Adjusted p-value threshold used when `method = "limma"`.
#' @param var_threshold Quantile threshold for MAD-based filtering when
#'   `method = "variance"` and `top_n = NULL`.
#' @param top_n Optional integer. Number of high-variance probes to retain per
#'   comparison when `method = "variance"`. If `NULL`, `var_threshold` is used.
#' @param variance_selection_mode Optional string recording whether variance
#'   filtering used `"top_n"` or `"threshold"` mode.
#' @param save_path Optional path to save the full result object as an `.rds` file.
#'
#' @return A nested list containing:
#'   \describe{
#'     \item{group/comparison entries}{Per-comparison filtered outputs, including
#'       `filtered_matrix`, `filtered_probes`, optional limma tables, and local metadata.}
#'     \item{`__summary__`}{Cross-comparison probe reuse summaries.}
#'     \item{`__metadata__`}{Global provenance metadata for the filtering run.}
#'   }
run_grouped_probe_filtering <- function(matrix_list,
                                        comparison_map,
                                        chip_id,
                                        method = filter_method,
                                        logfc_cutoff = logfc_cutoff,
                                        pval_cutoff = pval_cutoff,
                                        var_threshold = var_threshold,
                                        top_n = NULL,
                                        variance_selection_mode = NULL,
                                        save_path = NULL) {
  results <- list()
  probe_tracker <- list()
  
  for (group in names(comparison_map)) {
    results[[group]] <- list()
    
    for (label in names(comparison_map[[group]])) {
      pair <- comparison_map[[group]][[label]]
      ctrl <- paste0("m_", make.names(pair[1]))
      case <- paste0("m_", make.names(pair[2]))
      
      if (!(ctrl %in% names(matrix_list)) ||
          !(case %in% names(matrix_list))) {
        warning(glue::glue("Skipping comparison '{label}' — matrix not found."))
        next
      }
      
      ctrl_mat <- matrix_list[[ctrl]]
      case_mat <- matrix_list[[case]]
      combined <- cbind(ctrl_mat, case_mat)
      
      group_labels <- factor(c(rep("control", ncol(ctrl_mat)), rep("case", ncol(case_mat))))
      
      if (method == "limma") {
        # ---- Fit limma model with simple design ----
        design <- model.matrix( ~ group_labels)
        fit <- limma::lmFit(combined, design)
        fit <- limma::eBayes(fit)
        
        # ---- Extract top table and add probe_id column ----
        tab <- limma::topTable(fit,
                               coef = 2,
                               number = Inf,
                               sort.by = "none") |>
          tibble::rownames_to_column(var = "probe_id")
        
        # ---- Apply filtering ----
        filtered <- tab |>
          dplyr::filter(abs(logFC) >= logfc_cutoff, adj.P.Val < pval_cutoff)
        
        probe_ids <- filtered$probe_id
        probe_tracker[[paste(group, label, sep = "::")]] <- probe_ids
        
        results[[group]][[label]] <- list(
          filtered_probes = probe_ids,
          filtered_matrix = combined[probe_ids, , drop = FALSE],
          stats_table = filtered,
          limma_full = tab,
          metadata = list(
            chip = chip_id,
            group = group,
            label = label,
            comparison = pair
          )
        )
        
      } else if (method == "variance") {
        probe_ids <- select_high_variance_probes(
          combined,
          method = "mad",
          threshold = var_threshold,
          top_n = top_n
        )
        probe_tracker[[paste(group, label, sep = "::")]] <- probe_ids
        
        results[[group]][[label]] <- list(
          filtered_probes = probe_ids,
          filtered_matrix = combined[probe_ids, , drop = FALSE],
          metadata = list(
            chip = chip_id,
            group = group,
            label = label,
            comparison = pair
          )
        )
        
      } else {
        stop("Unknown filtering method: use 'limma' or 'variance'.")
      }
    }
  }
  
  # ------------------------------------------------------------------------------
  # SUMMARY: Cross-comparison probe usage
  # ------------------------------------------------------------------------------
  
  all_probes <- unlist(probe_tracker)
  
  results$`__summary__` <- list(
    probe_hits_by_group = probe_tracker,
    probe_multi_hit_counts = sort(table(all_probes), decreasing = TRUE)
  )
  
  # ------------------------------------------------------------------------------
  # METADATA: Pipeline provenance (global)
  # ------------------------------------------------------------------------------
  
  results$`__metadata__` <- list(
    pipeline_stage = "comparison_aware_probe_filtering",
    wrapper = "R/wrappers/run_grouped_probe_filtering.R",
    chip_id = chip_id,
    created_at = as.character(Sys.time()),
    method = method,
    variance_selection_mode = variance_selection_mode,
    parameters = list(
      logfc_cutoff = if (identical(method, "limma"))
        logfc_cutoff
      else
        NULL,
      pval_cutoff  = if (identical(method, "limma"))
        pval_cutoff
      else
        NULL,
      var_threshold = if (identical(method, "variance") &&
                          is.null(top_n))
        var_threshold
      else
        NULL,
      top_n = if (identical(method, "variance"))
        top_n
      else
        NULL
    )
  )
  
  # ------------------------------------------------------------------------------
  # NOTE (Design / Future Refactor)
  #
  # The current results object mixes three distinct concerns in a single structure:
  #   (1) Biological results (group → comparison → filtered outputs)
  #   (2) Cross-comparison summaries (__summary__)
  #   (3) Pipeline provenance metadata (__metadata__)
  #
  # This is functionally sufficient but structurally impure.
  #
  # Future refactor should separate these into:
  #
  #   list(
  #     metadata = ...,   # full provenance of filtering method and parameters
  #     data     = ...,   # nested biological results (current structure)
  #     summary  = ...    # cross-comparison probe reuse statistics
  #   )
  #
  # This will improve:
  #   - clarity of downstream access patterns
  #   - robustness of reporting pipelines (Quarto, etc.)
  #   - reproducibility and auditability
  #
  # IMPORTANT: This refactor will break downstream code and should only be done
  # in coordination with updates to analysis and reporting pipelines.
  # ------------------------------------------------------------------------------
  
  # ------------------------------------------------------------------------------
  # SAVE: Persist results object
  # ------------------------------------------------------------------------------
  
  if (!is.null(save_path)) {
    saveRDS(results, file = save_path)
  }
  
  return(results)
}
