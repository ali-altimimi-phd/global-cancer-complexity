# ------------------------------------------------------------------------------
# File: 02_run_analysis_pipeline.R
# Purpose: Execute the analysis pipeline
# Role: Top-level analysis runner
# Pipeline: Analysis
# Project: Cancer Complexity Analysis
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------
NULL
#' Run Analysis Pipeline for Global Cancer Microarray Data
#'
#'
#' This script is intended to be run **after** the global-cancer-complexity preprocessing pipeline.

run_analysis_pipeline <- function() {
  # ---- Stage 1: Configuration + Logging ----
  source(here::here("R/config/global_cancer/analysis_config.R"), local = FALSE)
  
  source(here::here("R/helpers/pipeline_logger.R"))
  logger <- start_log(analysis_pipeline_logfile)
  logger$log("🚀 Starting Global Cancer analysis pipeline...")
  
  # ---- Stage 2: Build matrix + comparison maps ----
  source(here::here("R/helpers/load_annotations.R"))
  load_annotations_if_needed(logger = logger)
  
  if (build_matrix_maps) {
    logger$log("⚙️ Building expression matrices and comparison maps...")
    
    # Use existing eset_list if available, otherwise load from disk
    if (exists("eset_list", inherits = TRUE)) {
      logger$log("⏭️ Using existing 'eset_list' from memory.")
    } else {
      if (!file.exists(eset_path))
        stop(glue::glue("❌ Missing ExpressionSets: {eset_path}"))
      
      loaded_objs <- load(eset_path)
      if (!"eset_list" %in% loaded_objs)
        stop("❌ 'eset_list' not found in loaded RData file.")
      
      eset_list <- get("eset_list")
      logger$log("✅ ExpressionSets loaded from disk.")
    }
    
    source(here::here("R/helpers/matrix_builders.R"), local = TRUE)
    
    matrices_hu35ksuba <- build_matrix_lists_by_tissue(eset_list$hu35ksuba)
    matrices_hu6800    <- build_matrix_lists_by_tissue(eset_list$hu6800)
    
    comparison_map_hu35ksuba <- define_predefined_comparisons(names(matrices_hu35ksuba))
    comparison_map_hu6800    <- define_predefined_comparisons(names(matrices_hu6800))
    
    dir.create(dirname(matrices_path), recursive = TRUE, showWarnings = FALSE)
    
    save(
      matrices_hu35ksuba,
      matrices_hu6800,
      comparison_map_hu35ksuba,
      comparison_map_hu6800,
      file = matrices_path
    )
    
    logger$log("💾 Matrix and comparison maps saved.")
  }
  
  # ---- Stage 3: Run group-aware filtering ----
  if (run_filtering) {
    logger$log("🔀 Checking for expression matrices and comparison maps...")
    source(here::here("R/helpers/load_results_into_global.R"))
    load_matrix_maps(matrices_path, overwrite = TRUE)
    
    logger$log("🧬 Running comparison-aware probe selection...")
    
    source(here::here("R/wrappers/run_grouped_probe_filtering.R"))
    source(here::here("R/filters/select_high_variance_probes.R"))
    
    dir.create(filtered_probes_dir, recursive = TRUE, showWarnings = FALSE)
    
    # ---- Validate variance selection mode ----
    if (!variance_selection_mode %in% c("top_n", "threshold")) {
      logger$log(
        sprintf(
          "⚠️ Invalid variance_selection_mode ('%s'); defaulting to 'threshold'.",
          variance_selection_mode
        ),
        section = "WARNING"
      )
      variance_selection_mode <- "threshold"
    }
    
    res_hu35ksuba <<- run_grouped_probe_filtering(
      matrix_list = matrices_hu35ksuba,
      comparison_map = comparison_map_hu35ksuba,
      chip_id = "hu35ksuba",
      method = filter_method,
      logfc_cutoff = logfc_cutoff,
      pval_cutoff = pval_cutoff,
      var_threshold = var_threshold,
      save_path = here::here(filtered_probes_dir, "filtered_probes_hu35ksuba.rds")
    )
    
    res_hu6800 <<- run_grouped_probe_filtering(
      matrix_list = matrices_hu6800,
      comparison_map = comparison_map_hu6800,
      chip_id = "hu6800",
      method = filter_method,
      logfc_cutoff = logfc_cutoff,
      pval_cutoff = pval_cutoff,
      var_threshold = var_threshold,
      save_path = here::here(filtered_probes_dir, "filtered_probes_hu6800.rds")
    )
    
    logger$log("✅ Completed probe selection.")
  }
  
  # ---- Stage 4: Pairwise comparisons ----
  source(here::here("R/helpers/load_annotations.R"))
  
  if (isTRUE(run_pairwise)) {
    source(here::here("R/helpers/load_results_into_global.R"))
    
    logger$log("🔀 Checking for comparison maps...", section = "PAIRWISE")
    load_matrix_maps(matrices_path, overwrite = TRUE)
    
    logger$log("🔀 Checking for filtered probes...", section = "PAIRWISE")
    load_filtered_results(filtered_probes_dir, overwrite = TRUE)
    
    load_annotations_if_needed(logger = logger)
    
    logger$log("🔀 Starting pairwise comparisons...", section = "PAIRWISE")
    source(here::here("R/wrappers/run_all_pairwise_comparisons.R"))
    
    for (gene_mode in gene_set_modes) {
      logger$log(
        sprintf("🔀 Running pairwise comparisons for gene_set_mode = %s", gene_mode),
        section = "PAIRWISE"
      )
      
      run_all_pairwise_comparisons(
        chips = chips,
        annotations = annotations,
        engines = engines,
        gene_set_mode = gene_mode,
        quantile_cutoff = quantile_cutoff,
        min_probes = min_probes
      )
    }
  }
  
  # ---- Stage 5: Aggregation + gene set annotation ----
  if (run_aggregator) {
    logger$log("🔀 Checking for pairwise comparison results...")
    source(here::here("R/helpers/load_results_into_global.R"))
    load_comparison_results(data_dir, overwrite = TRUE)
    
    logger$log("📦 Ensuring annotations are loaded...")
    source(here::here("R/helpers/load_annotations.R"))
    load_annotations_if_needed(logger = logger)
    
    logger$log("🔀 Starting results aggregation...")
    source(here::here("R/aggregate/aggregate_engine_results_by_engine.R"))
    source(here::here("R/helpers/gene_set_tools.R"))
    source(here::here("R/helpers/clean_aggregated_results.R"))
    
    agg_complexity <<- aggregate_engine_results_by_engine("complexity")
    agg_entropy    <<- aggregate_engine_results_by_engine("entropy")
    
    agg_complexity <<- attach_gene_set_names(agg_complexity, annotations)
    agg_entropy    <<- attach_gene_set_names(agg_entropy, annotations)
    
    complexity_df <<- clean_aggregated_results(agg_complexity)
    entropy_df    <<- clean_aggregated_results(agg_entropy)
    
    readr::write_csv(complexity_df,
                     file.path(aggregate_dir, "complexity_aggregated_results.csv"))
    readr::write_csv(entropy_df,
                     file.path(aggregate_dir, "entropy_aggregated_results.csv"))
    
    saveRDS(complexity_df,
            file.path(aggregate_dir, "complexity_aggregated_results.rds"))
    saveRDS(entropy_df,
            file.path(aggregate_dir, "entropy_aggregated_results.rds"))
    
    logger$log("✅ Aggregated complexity/entropy results saved.")
  }
  
  # ---- Stage 6: Summarize comparison-level patterns ----
  if (run_comparison_summary) {
    logger$log("🔀 Checking for cleaned pairwise comparison results...")
    source(here::here("R/helpers/load_results_into_global.R"))
    load_cleaned_results(aggregate_dir, overwrite = TRUE)
    
    logger$log("📊 Summarizing entropy + complexity across comparisons...")
    
    source("R/summarize/summarize_complexity.R")
    source("R/summarize/summarize_entropy.R")
    source("R/summarize/combine_entropy_complexity_summaries.R")
    
    summary_complexity_df <<- recalculate_complexity_summary(complexity_df)
    summary_entropy_df    <<- recalculate_entropy_summary(entropy_df)
    
    summaries_combined_df <<- combine_entropy_complexity_summaries(
      summary_complexity_df, summary_entropy_df
    )
    
    saveRDS(summaries_combined_df,
            file.path(summaries_dir, "summaries_combined_df.rds"))
    readr::write_csv(summaries_combined_df,
                     file.path(summaries_dir, "summaries_combined_df.csv"))
    
    logger$log("✅ Comparison-level summary table saved.")
  }
  

  # ---- Stage 7: Latent-space Python + notebook execution ----
  if (isTRUE(run_latent_python_scripts) || isTRUE(run_latent_notebooks)) {
    logger$log("🐍 Running latent-space analysis stages...", section = "LATENT")
    source(here::here("R/wrappers/run_latent_python_scripts.R"))
    source(here::here("R/wrappers/run_latent_notebooks.R"))
    
    if (isTRUE(run_latent_python_scripts)) {
      run_latent_python_scripts_stage(logger = logger)
    }
    
    if (isTRUE(run_latent_notebooks)) {
      run_latent_notebooks_stage(logger = logger)
    }
    
    logger$log("✅ Latent-space stages completed successfully.", section = "LATENT")
  } else {
    logger$log("⏭️ Latent-space stages skipped.", section = "LATENT")
  }

  # ---- Stage 8: Latent-space postprocessing ----
  if (isTRUE(run_latent_postprocessing)) {
    logger$timed("Latent-space postprocessing", {
      logger$log("🧾 Running latent-space postprocessing...", section = "LATENT_POST")
      
      source(here::here("R", "latent-space", "build_latent_postprocessing_outputs.R"), local = TRUE)
      
      run_latent_postprocessing(
        logger = logger,
        latent_output_dir = latent_output_dir,
        aggregate_dir = aggregate_dir,
        table_dir = latent_postprocessing_table_dir,
        plot_dir = latent_postprocessing_plot_dir,
        chip_id = latent_postprocessing_chip_id,
        mode = latent_postprocessing_mode,
        gene_set_name = latent_postprocessing_gene_set_name,
        complexity_delta_col = latent_postprocessing_complexity_delta_col
      )
      
      logger$log("✅ Latent-space postprocessing completed successfully.", section = "LATENT_POST")
    })
  } else {
    logger$log("⏭️ Latent-space postprocessing skipped.", section = "LATENT_POST")
  }

  # ---- Final message ----
  logger$log("🎉 Analysis pipeline completed successfully.", section = "PIPELINE")
}

if (sys.nframe() == 0) {
  run_analysis_pipeline()
}
