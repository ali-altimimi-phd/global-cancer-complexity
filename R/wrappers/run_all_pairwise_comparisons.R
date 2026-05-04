# ------------------------------------------------------------------------------
# File: run_all_reports.R
# Purpose: Generate comparison-level Quarto source reports and supporting
#   reporting artifacts by iterating over biological contrasts and invoking
#   data preparation, plot generation, and table export helpers.
# Role: Pipeline coordinator (batch execution wrapper)
# Pipeline: Reporting
# Project: Global Cancer Complexity
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------

#' Generate All Comparison-Level Quarto Report Sources
#'
#' This function prepares RDS data files, regenerates plots, writes
#' parameterized `.qmd` Quarto files, and delegates table writing to a helper
#' for each unique comparison in the data.
#'
#' Each generated `.qmd` file is placed in:
#' `quarto/reports/generated_reports/`, and can be included in the Quarto book.
#'
#' Plots are written to:
#' `quarto/resources/plots/`, organized by type (e.g., `probes/`, `go/`, `kegg/`).
#'
#' Reporting tables are generated via `write_all_comparison_tables()`, located in:
#' `R/helpers/write_all_comparison_tables.R`.
#'
#' @param summaries_combined_df Data frame containing summary results.
#' @param complexity_df Data frame of complexity statistics.
#' @param entropy_df Data frame of entropy statistics.
#' @param res_hu35ksuba Results object for hu35ksuba chip.
#' @param res_hu6800 Results object for hu6800 chip.
#' @param template_path Path to the Quarto template `.qmd` file.
#' @param rds_output_dir Directory to save RDS files.
#' @param plot_output_dir Directory to save plots.
#' @param qmd_output_dir Directory to save generated `.qmd` files.
#'
#' @return Invisibly returns \code{NULL}; writes report artifacts to disk.
run_all_pairwise_comparisons <- function(chips,
                                         annotations,
                                         engines,
                                         gene_set_mode,
                                         quantile_cutoff,
                                         min_probes) {
  # ---- Load dependencies ----
  source(here::here("R/helpers/pipeline_logger.R"))
  source(here::here("R/helpers/build_comparison_input_list.R"))
  source(here::here("R/helpers/gene_set_tools.R"))
  source(here::here("R/engines/complexity/compare_pair_complexity.R"))
  source(here::here("R/engines/entropy/compare_pair_entropy.R"))
  source(here::here("R/wrappers/run_pairwise_analysis.R"))

  # ---- Validate inputs ----
  valid_modes <- c("FULL", "GO_BP", "GO_MF", "KEGG", "MSIGDB")
  normalized_mode <- toupper(gene_set_mode)

  if (!normalized_mode %in% valid_modes) {
    stop(
      sprintf(
        "gene_set_mode must be one of %s. Received: %s",
        paste(valid_modes, collapse = ", "),
        gene_set_mode
      ),
      call. = FALSE
    )
  }

  if (!all(engines %in% c("complexity", "entropy"))) {
    stop("engines must contain only 'complexity' and/or 'entropy'.", call. = FALSE)
  }

  # ---- Logging setup ----
  pairwise_log_dir <- here::here(logs_dir, "pairwise")
  dir.create(pairwise_log_dir, recursive = TRUE, showWarnings = FALSE)

  logfile_path <- file.path(
    pairwise_log_dir,
    glue::glue("run_log_{tolower(normalized_mode)}_{format(Sys.time(), '%Y%m%d_%H%M%S')}.txt")
  )

  logger <- start_log(logfile = logfile_path)
  logger$log(
    glue::glue("🔀 Running pairwise comparisons for gene_set_mode = {normalized_mode}"),
    section = "PAIRWISE"
  )

  # ---- Ensure result directory exists ----
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- Precompute gene-set probe indices once per chip ----
  logger$log("⚙️ Building gene-set probe indices...", section = "PAIRWISE")
  gene_set_indices <- lapply(annotations, build_gene_set_probe_index)
  logger$log("✅ Gene-set probe indices built.", section = "PAIRWISE")

  # ---- Resolve gene-set source and ontology ----
  if (normalized_mode %in% c("GO_BP", "GO_MF")) {
    ontology <- sub("GO_", "", normalized_mode)
    source_key <- "GO"
  } else {
    ontology <- NULL
    source_key <- normalized_mode
  }

  # ---- Main loop: engine × chip ----
  for (engine in engines) {
    for (chip in chips) {
      logger$log(
        glue::glue("🔄 Starting {engine} comparison for chip: {chip} | Mode: {normalized_mode}"),
        section = "PAIRWISE"
      )

      filtered_name <- glue::glue("res_{chip}")
      comparison_name <- glue::glue("comparison_map_{chip}")

      if (!exists(filtered_name, envir = .GlobalEnv)) {
        stop(sprintf("Missing filtered results object: %s", filtered_name), call. = FALSE)
      }

      if (!exists(comparison_name, envir = .GlobalEnv)) {
        stop(sprintf("Missing comparison map object: %s", comparison_name), call. = FALSE)
      }

      if (!chip %in% names(annotations)) {
        stop(sprintf("Missing annotations for chip: %s", chip), call. = FALSE)
      }

      filtered_results <- get(filtered_name, envir = .GlobalEnv)
      comparison_map <- get(comparison_name, envir = .GlobalEnv)
      annotation_set <- annotations[[chip]]

      comparison_list <- build_comparison_input_list(
        comparison_map = comparison_map,
        filtered_results = filtered_results,
        chip_id = chip
      )

      # FULL means use the entire filtered comparison-specific probe space.
      # Gene-set modes use filtered probes intersected with eligible gene-set probes.
      gene_sets <- if (identical(source_key, "FULL")) {
        "FULL"
      } else {
        adaptive_gene_set_filter(
          annotation = annotation_set,
          source = source_key,
          quantile_cutoff = quantile_cutoff,
          min_probes = min_probes,
          ontology = ontology
        )
      }

      result <- logger$timed(glue::glue("{engine} - {chip} - {normalized_mode}"), {
        run_pairwise_analysis(
          comparison_list = comparison_list,
          gene_sets = gene_sets,
          engine = engine,
          gene_set_indices = gene_set_indices,
          verbose = FALSE
        )
      })

      suffix <- tolower(gsub("[^a-zA-Z0-9]+", "_", normalized_mode))
      out_name <- glue::glue("{engine}_results_{chip}_{suffix}")
      out_path <- here::here(data_dir, glue::glue("{out_name}.rds"))

      saveRDS(result, out_path)
      assign(out_name, result, envir = .GlobalEnv)

      logger$log(
        glue::glue("💾 Saved {engine} results for {chip} → {basename(out_path)}"),
        section = "PAIRWISE"
      )
    }
  }

  logger$log("✅ Pairwise comparisons completed.", section = "PAIRWISE")
  invisible(TRUE)
}
