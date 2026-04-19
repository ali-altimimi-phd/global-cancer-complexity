#' This script performs comparison-specific reduction of expression space 
#' by identifying probes most relevant to each tissue-specific normal-to-tumor 
#' transition. 
#' 
#' In doing so, it defines the subset of molecular features over which 
#' downstream complexity and entropy are evaluated, 
#' thereby focusing analysis on biologically salient transformation signals
#' rather than the full undifferentiated probe universe.
#' 
#' The resulting filtered probe sets preserve the local structure of each oncogenic
#' trajectory while also enabling comparison of recurrent probe usage across cancer types.
NULL
#' Run group-aware probe filtering by comparison and chip
#'
#' Executes probe filtering for all predefined group comparisons on a given chip.
#' Supports limma-based differential expression filtering or MAD/SD-based variability filtering.
#' Preserves per-group metadata and provides a summary of probe reuse across comparisons.
#'
#' @param matrix_list A named list of expression matrices from `build_matrix_lists_by_tissue()`
#' @param comparison_map A named list of group → label → tissue comparisons
#' @param chip_id A string identifying the chip (e.g., "hu35ksuba")
#' @param method Filtering method: either `"limma"` or `"variance"`
#' @param logfc_cutoff Log2 fold-change threshold (for limma only)
#' @param pval_cutoff Adjusted p-value threshold (for limma only)
#' @param var_threshold Quantile threshold for MAD/SD filtering (for variance only)
#' @param save_path Optional. Path to save the full result object as `.rds`
#'
#' @return A nested list structured as results[group][[label]], each containing:
#'   - `filtered_matrix`: matrix of selected probes
#'   - `filtered_probes`: vector of probe IDs
#'   - `stats_table`: differential expression results (limma only)
#'   - `limma_full`: full limma table (limma only)
#'   - `metadata`: chip/group/label information
#'   - `__summary__`: list of all probes per comparison and counts across groups
run_grouped_probe_filtering <- function(matrix_list,
                                        comparison_map,
                                        chip_id,
                                        method = filter_method,
                                        logfc_cutoff = logfc_cutoff,
                                        pval_cutoff = pval_cutoff,
                                        var_threshold = var_threshold,
                                        save_path = NULL) {
  results <- list()
  probe_tracker <- list()
  
  for (group in names(comparison_map)) {
    results[[group]] <- list()
    
    for (label in names(comparison_map[[group]])) {
      pair <- comparison_map[[group]][[label]]
      ctrl <- paste0("m_", make.names(pair[1]))
      case <- paste0("m_", make.names(pair[2]))
      
      if (!(ctrl %in% names(matrix_list)) || !(case %in% names(matrix_list))) {
        warning(glue::glue("Skipping comparison '{label}' — matrix not found."))
        next
      }
      
      ctrl_mat <- matrix_list[[ctrl]]
      case_mat <- matrix_list[[case]]
      combined <- cbind(ctrl_mat, case_mat)
      
      group_labels <- factor(c(
        rep("control", ncol(ctrl_mat)),
        rep("case",    ncol(case_mat))
      ))
      
      if (method == "limma") {
        # ---- Fit limma model with simple design ----
        design <- model.matrix(~ group_labels)
        fit <- limma::lmFit(combined, design)
        fit <- limma::eBayes(fit)
        
        # ---- Extract top table and add probe_id column ----
        tab <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none") |>
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
        probe_ids <- select_high_variance_probes(combined,
                                                 method = "mad",
                                                 threshold = var_threshold,
                                                 top_n = top_n)
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
  
  # ---- Add summary of probe hits across comparisons ----
  all_probes <- unlist(probe_tracker)
  results$`__summary__` <- list(
    probe_hits_by_group = probe_tracker,
    probe_multi_hit_counts = sort(table(all_probes), decreasing = TRUE)
  )
  
  # ---- Optionally save to disk ----
  if (!is.null(save_path)) {
    saveRDS(results, file = save_path)
  }
  
  return(results)
}
