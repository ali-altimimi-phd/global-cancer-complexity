# ------------------------------------------------------------------------------
# File: preprocessing_config.R
# Purpose: Configuration for the global cancer preprocessing pipeline
# Role: Defines toggles, download settings, metadata corrections, chip mappings,
#       file paths, logging destinations, and ML/latent-space export settings.
# Pipeline: Preprocessing
# Project: Cancer Complexity Analysis
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------

# ==============================================================================
# 1. Study identifier
# ==============================================================================

# Used to define study-specific input/output directory structure.
study_name <- "global_cancer"


# ==============================================================================
# 2. Main preprocessing pipeline toggles
# ==============================================================================

# Controls major preprocessing stages.
download_enabled <- FALSE
build_esets      <- FALSE
process_metadata <- FALSE
run_annotation   <- FALSE


# ==============================================================================
# 3. ML / latent-space export toggles and parameters
# ==============================================================================

# Exports ExpressionSet-derived files for downstream Python/Jupyter workflows.
# This stage does not train the VAE and does not perform latent-space analysis.
export_ml_inputs <- TRUE
run_ml_pca_check <- TRUE

# Current latent-space work uses hu35ksuba.
# Later this can become names(geo_chip_map) or a vector such as:
#   c("hu35ksuba", "hu6800")
ml_chip_id       <- "hu35ksuba"
ml_filter_method <- "variance"
ml_top_n         <- 3000


# ==============================================================================
# 4. Download settings
# ==============================================================================

# If TRUE, report download targets without downloading files.
# Use FALSE for actual downloading.
dry_run <- TRUE

# Limit number of files downloaded during testing.
# Use NULL to download all files.
download_limit <- 3

# Remote source for CEL files.
ftp_base <- "ftp://ftp.ebi.ac.uk/biostudies/fire/E-GEOD-/928/E-GEOD-68928/Files/"

# GEO accession used for metadata retrieval.
geo_accession <- "GSE68928"


# ==============================================================================
# 5. Chip/platform mappings
# ==============================================================================

# Maps internal chip identifiers to GEO platform identifiers.
geo_chip_map <- list(
  hu35ksuba = "GPL98",
  hu6800    = "GPL80"
)


# ==============================================================================
# 6. Metadata normalization fixes
# ==============================================================================

# Dataset-specific tissue-label corrections.
fix_tissue_labels <- c(
  "organism part: Kideny"         = "organism part: Kidney",
  "organism part: Lymphod Tissue" = "organism part: Lymphoid Tissue",
  "organism part: Lymphoid"       = "organism part: Lymphoid Tissue"
)

# Dataset-specific disease-label corrections.
fix_disease_labels <- c(
  "disease state: large-Bcell lymphoma"               = "disease state: large B-cell lymphoma",
  "disease state: bladder transitonal cell carcinoma" = "disease state: bladder transitional cell carcinoma"
)


# ==============================================================================
# 7. Core input, output, and logging paths
# ==============================================================================

# Raw/local input directories.
local_cel_dir <- here::here("data", study_name, "CEL")
local_geo_dir <- here::here("data", study_name, "GEO", geo_accession)

# Preprocessing logs and RData output directories.
logs_dir <- here::here("output", study_name, "logs", "preprocess")
data_dir <- here::here("output", study_name, "RData")

# Avoids accidental double nesting such as:
#   output/global_cancer/logs/preprocess/preprocess/...
preprocess_pipeline_logfile <- file.path(
  logs_dir,
  "preprocess_pipeline_log.txt"
)

# Core preprocessing outputs.
eset_path <- file.path(
  data_dir,
  paste0(study_name, "_eset_list.RData")
)

annotations_path <- file.path(
  data_dir,
  "annotations",
  "full_chip_annotations.rds"
)


# ==============================================================================
# 8. ML / latent-space export output paths
# ==============================================================================

# Canonical downstream inputs for Python scripts, notebooks, and ML workflows.
ml_output_dir <- here::here(
  "data",
  study_name,
  "processed",
  "ml_inputs"
)

# PCA sanity-check and latent-related preprocessing plots.
# Currently routed to Quarto resources so they can be included directly in reports.
ml_plot_dir <- here::here(
  "quarto",
  "resources",
  "plots",
  "latent"
)


# ==============================================================================
# 9. Ensure required directories exist
# ==============================================================================

dir.create(dirname(preprocess_pipeline_logfile), recursive = TRUE, showWarnings = FALSE)

dir.create(local_cel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(local_geo_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(dirname(eset_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(annotations_path), recursive = TRUE, showWarnings = FALSE)

dir.create(ml_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ml_plot_dir, recursive = TRUE, showWarnings = FALSE)