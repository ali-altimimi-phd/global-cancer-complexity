#' This script operationalizes pairwise biological inference across 
#' tissue-specific normal-to-malignant transitions by applying complexity 
#' and entropy engines to chip-specific filtered comparison inputs. 
#' 
#' It supports both whole-probe and gene set–restricted analyses, 
#' with adaptive gene-set eligibility determined by platform-specific 
#' annotation coverage.
#'
#' In doing so, it constitutes the core execution layer through which malignant
#' transformation is assessed both as a change in overall transcriptomic 
#' complexity and as a redistribution of complexity across functional 
#' biological processes.
NULL

# R/wrappers/run_all_pairwise_comparisons.R

#' Run Pairwise Comparisons for All Chips and Engines
#'
#' Executes complexity- or entropy-based comparisons across all chip platforms,
#' using group-specific filtered matrices and adaptively filtered gene sets.
#'
#' Gene sets (GO, KEGG, MSigDB) are selected adaptively per chip using
#' platform-specific annotation coverage and probe-count thresholds.
#' Results are saved as `.rds` objects and summarized via structured logs.
#'
#' @param chips Character vector of chip IDs (e.g., `c("hu35ksuba", "hu6800")`)
#' @param annotations List of gene set annotations by chip
#' @param engines Character vector: `"complexity"` and/or `"entropy"`
#' @param gene_set_mode One of `"ALL"`, `"GO_BP"`, `"GO_MF"`, `"KEGG"`, `"MSIGDB"`
#' @param quantile_cutoff Numeric (e.g., 0.75). Quantile cutoff for adaptive filtering
#' @param min_probes Integer. Minimum number of probes per gene set
#'
#' @return No return value. Saves `.rds` files and logs to disk.
#' @export
run_all_pairwise_comparisons <- function(
    chips,
    annotations,
    engines,
    gene_set_mode,
    quantile_cutoff,
    min_probes
) {
  # ---- Load dependencies ----
  source(here::here("R/helpers/pipeline_logger.R"))
  source(here::here("R/helpers/build_comparison_input_list.R"))
  source(here::here("R/helpers/gene_set_tools.R"))
  source(here::here("R/engines/complexity/compare_pair_complexity.R"))
  source(here::here("R/engines/entropy/compare_pair_entropy.R"))
  source(here::here("R/wrappers/run_pairwise_analysis.R"))
  
  # ---- Logging setup ----
  logfile_path <- here::here(logs_dir, "pairwise",
                             glue::glue("run_log_{format(Sys.time(), '%Y%m%d_%H%M%S')}.txt"))
  logger <- start_log(logfile = logfile_path)
  
  # ---- Main loop: chip × engine ----
  for (engine in engines) {
    for (chip in chips) {
      logger$log(glue::glue("🔄 Starting {engine} comparison for chip: {chip} | Mode: {gene_set_mode}"))
      
      # Lookup required inputs
      filtered_results <- get(glue::glue("res_{chip}"), envir = .GlobalEnv)
      comparison_map   <- get(glue::glue("comparison_map_{chip}"), envir = .GlobalEnv)
      annotation_set   <- annotations[[chip]]
      
      # Build comparison matrix pairs
      comparison_list <- build_comparison_input_list(
        comparison_map = comparison_map,
        filtered_results = filtered_results,
        chip_id = chip
      )
      
      # ---- Compute dynamic thresholds ----
      thresholds_by_chip <- compute_dynamic_gene_set_thresholds(
        annotations = annotations,
        output_log_file = here::here(logs_dir, "pairwise", "dynamic_thresholds.txt")
      )
      
      # ---- Adaptive gene set selection ----
      normalized_mode <- toupper(gene_set_mode)
      
      # Handle GO sub-ontologies
      if (normalized_mode %in% c("GO_BP", "GO_MF")) {
        ontology <- sub("GO_", "", normalized_mode)  # Extract "BP" or "MF"
        source_key <- "GO"
      } else {
        ontology <- NULL
        source_key <- normalized_mode
      }
      
      gene_sets <- if (source_key == "ALL") {
        "ALL"
      } else {
        adaptive_gene_set_filter(
          annotation = annotation_set,
          source = source_key,
          quantile_cutoff = thresholds_by_chip[[chip]][[source_key]]$quantile_cutoff,
          min_probes = thresholds_by_chip[[chip]][[source_key]]$min_probes,
          ontology = ontology
        )
      }
      
      # ---- Run complexity/entropy engine ----
      result <- logger$timed(glue::glue("{engine} - {chip} - {gene_set_mode}"), {
        run_pairwise_analysis(
          comparison_list = comparison_list,
          gene_sets = gene_sets,
          engine = engine,
          annotations = annotations,
          verbose = FALSE
        )
      })
      
      # ---- Save results ----
      suffix <- tolower(gsub("[^a-zA-Z0-9]+", "_", gene_set_mode))
      out_name <- glue::glue("{engine}_results_{chip}_{suffix}")
      out_path <- here::here(data_dir, glue::glue("{out_name}.rds"))
      
      saveRDS(result, out_path)
      assign(out_name, result, envir = .GlobalEnv)
      
      logger$log(glue::glue("💾 Saved {engine} results for {chip} → {basename(out_path)}"))
    }
  }
  
  logger$log("✅ All pairwise comparisons completed.")
}
