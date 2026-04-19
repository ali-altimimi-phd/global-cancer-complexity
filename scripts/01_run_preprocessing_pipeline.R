# ------------------------------------------------------------------------------
# File: 01_run_preprocessing_pipeline.R
# Purpose: Runs preprocessing pipeline
# Role: Runner
# Pipeline: N/A
# Project: Cancer Complexity Analysis
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------

#' Run Preprocessing Pipeline for Global Cancer Microarray Data
#'
#' Executes the preprocessing workflow for Affymetrix CEL files, including:
#' downloading, normalizing, cleaning metadata, assigning condition/tissue labels,
#' and annotating all chip-specific probe sets.
#'
#' This script is intended to be run **before** the analysis pipeline.
#' It writes cleaned `eset_list` (as .RData) and full-chip annotations (as .rds) to disk.
#'
#' @note The full pipeline may take up-to 75 minutes depending on file count and network speed.
#'
#' @return No return value. Writes output to disk.
#' @export
run_preprocessing_pipeline <- function() {
  # ---- Config ----
  source(here::here("R/config/global_cancer/preprocessing_config.R"))
  
  # ---- Logging ----
  source(here::here("R/helpers/pipeline_logger.R"))
  logger <- start_log(preprocess_pipeline_logfile)
  logger$log("🚀 Starting Preprocessing Pipeline...")
  
  # ---- Step 1: Download CEL files (optional) ----
  if (download_enabled) {
    logger$log("📥 Downloading CEL files via FTP...")
    source(here::here("R/preprocessing/download_cel_files.R"))
    downloaded <- download_cel_files(
      ftp_base = ftp_base,
      local_cel_dir = local_cel_dir,
      download_limit = download_limit,
      dry_run = dry_run,
      logger = logger
    )
    logger$log(glue::glue("✅ FTP download complete. {length(downloaded)} files processed."))
  } else {
    logger$log("⏭️ FTP download skipped (download_enabled = FALSE).")
  }

  # ---- Step 2: Build ExpressionSets from CEL files (optional) ----
  if (build_esets) {
    logger$log("🧬 Building ExpressionSet objects from CEL files...")
    source(here::here("R/preprocessing/build_expression_sets.R"))
    eset_list <- build_expression_sets(
      cel_dir = local_cel_dir,
      log_dir = logs_dir,
      logger = logger
    )
  } else {
    logger$log("⏭️ ExpressionSet building skipped (build_esets = FALSE).")
  }
  
  if (!build_esets) {
    load(eset_path)
  }
    
  # ---- Step 3: Attach and Process GEO Metadata ----
  if (process_metadata) {
    logger$log("🧬 Processing GEO metadata: attach → clean → label...")
    
    # Attach metadata
    source(here::here("R/preprocessing/attach_geo_metadata.R"))
    eset_list <- attach_geo_metadata(
      eset_list = eset_list,
      geo_dir = local_geo_dir,
      logger = logger
    )
    
    # Clean and label metadata
    source(here::here("R/helpers/clean_and_label_metadata.R"), local = TRUE)
    eset_list <- clean_and_label_metadata(
      eset_list = eset_list,
      tissue_fixes = fix_tissue_labels,
      disease_fixes = fix_disease_labels,
      logger = logger
    )
    
    logger$log("✅ GEO metadata processing complete.")
  } else {
    logger$log("⏭️ GEO metadata processing skipped (process_metadata = FALSE).")
  }
  
  # ---- Step 4: Annotate all chip probes ----
  if (run_annotation) {
    logger$log("🔎 Annotating all probes per chip...")
    source(here::here("R/wrappers/run_annotate_chip_probes.R"), local = TRUE)
    
    annotations <- run_annotate_chip_probes(
      eset_list        = eset_list,
      chip_ids         = names(geo_chip_map),
      annotations_path = annotations_path
    )
    
    logger$log("✅ Probe annotation complete.")
  } else {
    logger$log("⏭️ Annotation step skipped (run_annotation = FALSE).")
  }
  
  # ---- Step 5: Save ExpressionSets ----
  save(eset_list, file = eset_path)
  logger$log(glue::glue("💾 ExpressionSets saved to {basename(eset_path)}"))
  
  logger$log("🏁 Preprocessing pipeline completed successfully.")
}
