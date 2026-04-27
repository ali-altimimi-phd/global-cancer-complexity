# ------------------------------------------------------------------------------
# File: preprocessing_config.R
# Purpose: Define configuration parameters, file paths, toggles, chip mappings,
#   metadata sources, logging destinations, and ML/latent-space export settings
#   for the global cancer preprocessing pipeline
# Role: Preprocessing configuration file
# Pipeline: Preprocessing
# Project: Cancer Complexity Analysis
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------

# ---- Study Identifier ----
# Used to define directory structure for multiple data sources
study_name <- "global_cancer"

# ---- Preprocessing Pipeline Toggles ----
download_enabled      <- FALSE
build_esets           <- TRUE
process_metadata      <- TRUE
run_annotation        <- TRUE

# ---- ML / Latent-Space Export Toggles ----
# Exports ExpressionSet-derived files for downstream Python/Jupyter workflows.
# This stage does not train the VAE and does not perform latent-space analysis.
export_ml_inputs      <- TRUE
run_ml_pca_check      <- TRUE

# Current latent-space work uses hu35ksuba. Later this can become names(geo_chip_map)
# or a vector such as c("hu35ksuba", "hu6800").
ml_chip_id            <- "hu35ksuba"
ml_filter_method      <- "variance"
ml_top_n              <- 3000

# ---- Download Toggles ----
dry_run               <- TRUE      # Use false for actual downloading
download_limit        <- 3         # Use NULL for all files

# ---- FTP or GEO Source Info ----
ftp_base <- "ftp://ftp.ebi.ac.uk/biostudies/fire/E-GEOD-/928/E-GEOD-68928/Files/"

geo_accession <- "GSE68928"

geo_chip_map <- list(
  hu35ksuba = "GPL98",
  hu6800    = "GPL80"
)

# ---- Metadata Fixes ----
# These are specific to the global cancer data set
fix_tissue_labels <- c(
  "organism part: Kideny"         = "organism part: Kidney",
  "organism part: Lymphod Tissue" = "organism part: Lymphoid Tissue",
  "organism part: Lymphoid"       = "organism part: Lymphoid Tissue"
)
fix_disease_labels <- c(
  "disease state: large-Bcell lymphoma"               = "disease state: large B-cell lymphoma",
  "disease state: bladder transitonal cell carcinoma" = "disease state: bladder transitional cell carcinoma"
)

# ---- Local Paths ----
local_cel_dir <- here::here("data", study_name, "CEL")
local_geo_dir <- here::here("data", study_name, "GEO", "GSE68928")
logs_dir      <- here::here("output", study_name, "logs", "preprocess")
data_dir      <- here::here("output", study_name, "RData")

# ---- Output File Paths ----
# Avoids accidental double nesting: output/global_cancer/logs/preprocess/preprocess/...
preprocess_pipeline_logfile <- file.path(logs_dir, "preprocess_pipeline_log.txt")

eset_path        <- here::here(data_dir, paste0(study_name, "_eset_list.RData"))
annotations_path <- here::here(data_dir, "annotations", "full_chip_annotations.rds")

# ---- ML / Latent-Space Output Paths ----
# These are canonical downstream inputs for notebooks and other ML workflows.
ml_output_dir <- here::here("data", study_name, "processed", "ml_inputs")
# ml_plot_dir   <- here::here("output", study_name, "plots", "preprocessing")
ml_plot_dir   <- here::here("quarto", "resources", "plots", "latnt")

# ---- Ensure Directories Exist ----
dir.create(dirname(preprocess_pipeline_logfile), recursive = TRUE, showWarnings = FALSE)
dir.create(local_cel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(local_geo_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(eset_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(annotations_path), recursive = TRUE, showWarnings = FALSE)
dir.create(ml_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ml_plot_dir, recursive = TRUE, showWarnings = FALSE)
