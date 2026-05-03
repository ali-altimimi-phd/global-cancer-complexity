#' This script executes comparison-level inference across tissue-specific
#' normal-to-malignant transitions by applying the selected analytic engine
#' to either the full filtered probe space or to biologically defined gene-set
#' subspaces. It enforces minimum data sufficiency thresholds and records
#' skipped comparisons, thereby ensuring that downstream conclusions are drawn
#' only from statistically admissible matrix contrasts.
NULL
#' Run Pairwise Complexity or Entropy Comparisons
#'
#' Wrapper to execute pairwise comparisons using complexity or entropy engines.
#' Accepts group-aware filtered matrices and annotations. Skips comparisons with
#' insufficient probes and logs them for review.
#'
#' @param comparison_list List of comparison metadata and matrices (from build_comparison_input_list)
#' @param gene_sets Character vector of gene sets to test (e.g., "FULL", GO IDs, KEGG IDs)
#' @param engine Engine to use: "complexity" or "entropy"
#' @param annotations Named list of annotation tables by chip
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A named list of results, grouped by comparison and gene set.
#' @export
# run_pairwise_analysis <- function(
    #     comparison_list,
#     gene_sets = "FULL",
#     engine = c("complexity", "entropy"),
#     annotations = NULL,
#     verbose = TRUE
# ) {
run_pairwise_analysis <- function(comparison_list,
                                  gene_sets = "FULL",
                                  engine = c("complexity", "entropy"),
                                  gene_set_indices = NULL,
                                  verbose = TRUE) {
  engine <- match.arg(engine)
  
  if (!(is.character(gene_sets) && length(gene_sets) >= 1)) {
    stop("❌ gene_sets must be a character vector: 'FULL' or gene set IDs.")
  }
  
  results <- list()
  skipped <- list()
  
  for (cmp in comparison_list) {
    cmp_name <- cmp$name
    group    <- cmp$group
    chip     <- cmp$chip
    mat1     <- cmp$m1
    mat2     <- cmp$m2
    
    if (nrow(mat1) < 2 || nrow(mat2) < 2) {
      message(glue::glue(
        "⚠️ Skipping: {cmp_name} ({chip}) — <2 probes in one or both matrices"
      ))
      skipped[[cmp_name]] <- list(group = group,
                                  chip = chip,
                                  reason = "<2 probes")
      next
    }
    
    if (verbose)
      message(glue::glue("🔬 Running {engine} for: {cmp_name} ({chip})"))
    
    cmp_result <- list()
    
    for (gs in gene_sets) {
      id <- if (gs == "FULL")
        "FULL"
      else
        make.names(gs)
      
      # if (gs == "FULL") {
      #   probes <- NULL
      # } else {
      #   if (is.null(annotations) || is.null(annotations[[chip]])) {
      #     stop("❌ Annotations missing or chip not found: ", chip)
      #   }
      #   probes <- get_probes_for_set(gs, annotations[[chip]])
      # }
      
      if (gs == "FULL") {
        probes <- NULL
      } else {
        if (is.null(gene_set_indices) || is.null(gene_set_indices[[chip]])) {
          stop("❌ Gene-set index missing or chip not found: ", chip)
        }
        
        probes <- get_probes_for_set(gs, gene_set_indices[[chip]])
      }
      
      result <- switch(
        engine,
        "complexity" = compare_matrix_pair_complexity(
          label = cmp_name,
          chip = chip,
          mat1 = mat1,
          mat2 = mat2,
          filter_probes = probes
        ),
        "entropy" = compare_matrix_pair_entropy(
          label = cmp_name,
          chip = chip,
          mat1 = mat1,
          mat2 = mat2,
          filter_probes = probes
        ),
        stop("❌ Unknown engine: ", engine)
      )
      
      if (!is.null(result)) {
        cmp_result[[id]] <- result
        message(
          glue::glue(
            "✅ Completed: {cmp_name} | Chip: {chip} | Gene set: {gs} | Engine: {engine}"
          )
        )
      } else {
        message(
          glue::glue(
            "⚠️ Skipped: {cmp_name} | Chip: {chip} | Gene set: {gs} | Engine: {engine}"
          )
        )
      }
    }
    
    results[[cmp_name]] <- cmp_result
  }
  
  # ---- Save matrix-level skip log if needed ----
  if (length(skipped) > 0) {
    skip_df <- purrr::map_dfr(names(skipped), function(name) {
      tibble::tibble(
        comparison = name,
        group = skipped[[name]]$group,
        chip = skipped[[name]]$chip,
        reason = skipped[[name]]$reason
      )
    })
    
    log_dir <- here::here("output", "logs", "pairwise")
    if (!dir.exists(log_dir))
      dir.create(log_dir, recursive = TRUE)
    
    log_file <- file.path(log_dir, "skipped_comparisons.csv")
    readr::write_csv(skip_df, log_file)
    message(glue::glue("📝 Skipped comparison log saved to: {basename(log_file)}"))
  }
  
  return(results)
}
