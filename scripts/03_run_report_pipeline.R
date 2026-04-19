#' Run Global Cancer Expression Report Pipeline
#'
#' Report generation
#'
#' @note Requires datasets created via `run_analysis_pipeline()`.
#' @seealso \code{R/config/analysis_pipeline_config.R}
#' @return No return value; writes all outputs to disk.
#' @export
run_report_pipeline <- function() {
  # ---- Config ----
  source(here::here("R/config/global_cancer/report_config.R"),
         local = FALSE)
  
  # ---- Logging ----
  source(here::here("R/helpers/pipeline_logger.R"))
  logger <- start_log(reports_pipeline_logfile)
  logger$log("🚀 Starting Global Cancer report pipeline")
  
  # ---- Stage 1: Load data sets ----
  source(here::here("R/helpers/load_results_into_global.R"))
  
  logger$log("🔀 Checking for summary results...")
  load_summaries_df(summaries_dir, overwrite = TRUE)
  
  logger$log("🔀 Checking for cleaned pairwise comparison results...")
  load_cleaned_results(aggregate_dir, overwrite = TRUE)

  logger$log("🔀 Checking for filtered probes...")
  load_filtered_results(filtered_probes_dir, overwrite = TRUE)

  # ---- Stage 2: Create summary reports ----
  if (run_summary_reports) {
    logger$log("📄 Generating summary reports...")
    source(here::here("R/reports/summarize_pairwise_results.R"))
    
    entropy_report    <- summarize_pairwise_results(entropy_df, engine = "entropy")
    complexity_report <- summarize_pairwise_results(complexity_df, engine = "complexity")
    
    readr::write_lines(entropy_report,
                       file.path(reports_dir, "entropy_summary.txt"))
    readr::write_lines(complexity_report,
                       file.path(reports_dir, "complexity_summary.txt"))
    
    readr::write_csv(entropy_df,
                     file.path(reports_dir, "entropy_cleaned.csv"))
    readr::write_csv(complexity_df,
                     file.path(reports_dir, "complexity_cleaned.csv"))
    
    saveRDS(entropy_df,
            file.path(reports_dir, "entropy_cleaned.rds"))
    saveRDS(complexity_df,
            file.path(reports_dir, "complexity_cleaned.rds"))
    rmarkdown::render(here::here(("R/reports/pairwise_summary_report.Rmd")))
    
    logger$log("✅ Reports saved.")
    
  }
  
  # ---- Stage 3: Create reports from .qmd template ----
  
  # This might be a misnomer as we are not running reports but we are creating
  # comparison qmds from a template qmd
  if (run_reports) {
    source(here::here("R/wrappers/run_all_reports.R"))
    run_all_reports(
      summaries_combined_df = summaries_combined_df,
      complexity_df = complexity_df,
      entropy_df = entropy_df,
      res_hu35ksuba = res_hu35ksuba,
      res_hu6800 = res_hu6800
    )
  }
  
  # ---- Stage 4: GO clusters analysis ----
  if (run_go_clustering) {
    logger$log("🧬 Running GO clustering analysis...", {go_mode})
    
    source(here::here("R/wrappers/run_go_clustering.R")) # Dispatcher for all comparisons
    
    run_go_clustering_main(
      summaries_combined_df = summaries_combined_df,
      complexity_df = complexity_df,
      entropy_df = entropy_df,
      output_dir_base = here::here("quarto", "resources", "tables", "clusters"),
      go_mode = go_mode
    )
    
    logger$log("✅ GO clustering analysis complete.")
  }
  
  # ---- Stage 5: Create GO visualizations ----
  # source("R/visualizations/plot_entropy_go_dag.R")
  # plot_entropy_go_dag("PB/B-ALL")
  
  # ---- Stage 6: Render comparison reports from .qmd template ----
  if (run_quarto) {
    quarto::quarto_render("quarto")
  }

  # ---- Final message ----
  logger$log("🎉 Report pipeline completed successfully.")
}