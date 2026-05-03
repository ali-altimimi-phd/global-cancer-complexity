# ------------------------------------------------------------------------------
# File: analysis_config.R
# Purpose: Define configuration parameters, file paths, toggles, engines,
#   notebook execution settings, postprocessing stages, and logging destinations
#   for the global cancer analysis pipeline
# Role: Analysis configuration file
# Pipeline: Analysis
# Project: Cancer Complexity Analysis
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------

# ---- Study Identifier ----
# Used to define directory structure for multiple data sources
study_name <- "global_cancer"

# ---- Analysis Pipeline Toggles ----
build_matrix_maps <- FALSE
run_filtering <- FALSE
run_pairwise <- FALSE
run_aggregator <- TRUE
run_comparison_summary <- TRUE

# ---- Latent-Space Notebook Toggles ----
run_latent_notebooks <- FALSE

# If TRUE, notebooks are executed in place. This updates notebook outputs.
# If FALSE, executed copies are written to notebook_executed_dir.
execute_notebooks_inplace <- FALSE

# Set to -1 for no timeout; useful for VAE training notebooks.
notebook_timeout_seconds <- -1

latent_python_exe <- "C:/Users/drziy/anaconda3/envs/ml/python.exe"
jupyter_command <- latent_python_exe
jupyter_prefix_args <- c("-m", "jupyter")
jupyter_command <- latent_python_exe
jupyter_prefix_args <- c("-m", "jupyter")

latent_notebook_dir <- here::here("projects", "cancer-latent-space", "notebooks")
notebook_executed_dir <- here::here("output", study_name, "notebooks", "executed")

latent_notebooks <- c(
  "00_latent_space_config.ipynb",
  "01_load_and_inspect_data.ipynb",
  "02_prepare_data_for_vae.ipynb",
  "03_train_first_vae.ipynb",
  "04_latent_space_analysis_standardized.ipynb",
  "05_latent_complexity_analysis_standardized.ipynb",
  "06_between_cancer_geometry_in_latent_space_enhanced.ipynb"
)

# ---- Latent-Space Postprocessing Toggles ----
run_latent_postprocessing <- FALSE

latent_project_dir <- here::here("output", "global_cancer")
latent_postprocessing_table_dir <- file.path(latent_project_dir, "tables", 
                                             "latent")
latent_postprocessing_plot_dir <- file.path(latent_project_dir, "plots", 
                                            "latent")

latent_postprocessing_chip_id <- "hu35ksuba"
latent_postprocessing_mode <- "ALL"
latent_postprocessing_gene_set_name <- "ALL"

dir.create(latent_postprocessing_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(latent_postprocessing_plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Probe Selection Parameters ----

# Primary selection method:
#   "limma"    = differential-expression-based probe selection
#   "variance" = variability-based probe selection (see mode below)
filter_method <- "limma"

# ---- Limma mode parameters ----
# Used only when filter_method == "limma"
logfc_cutoff <- 0.33
pval_cutoff  <- 0.05

# ---- Variance mode parameters ----
# Used only when filter_method == "variance"

# Selection mode:
#   "top_n"     = retain a fixed number of most variable probes per comparison
#   "threshold" = retain probes at or above a variability quantile cutoff
#
# If an invalid value is supplied, the pipeline will log a warning and
# default to "threshold".
variance_selection_mode <- "top_n"

# Fixed-count selection (used only when variance_selection_mode == "top_n")
analysis_top_n <- 3000

# Threshold-based selection (used only when variance_selection_mode == "threshold")
# Must be between 0 and 1 (quantile of variability distribution)
var_threshold <- 0.75

# ---- Pairwise Comparison Settings ----
chips <- c("hu35ksuba", "hu6800")
geo_chip_map <- list(hu35ksuba = "GPL98", hu6800 = "GPL80")
engines <- c("complexity", "entropy")

# ---- Gene Set Filter ----
gene_set_modes <- c("FULL", "GO_BP", "GO_MF", "KEGG", "MSIGDB")

# For repair / partial reruns, use one or more modes, e.g.:
# gene_set_modes <- c("KEGG")
# gene_set_modes <- c("GO_BP", "GO_MF")

quantile_cutoff <- 0.75
min_probes <- 5

# ---- Local Paths ----
logs_dir <- here::here("output", study_name, "logs")
analysis_logs_dir <- file.path(logs_dir, "analysis")

data_dir <- here::here("output", study_name, "RData")
filtered_probes_dir <- here::here(data_dir, "filtered_probes")
aggregate_dir <- here::here(data_dir, "aggregated")
summaries_dir <- here::here(data_dir, "summaries")
reports_dir <- here::here("output", study_name, "reports")

# ---- Output File Paths ----
analysis_pipeline_logfile <- file.path(analysis_logs_dir, "analysis_pipeline_log.txt")
eset_path <- here::here(data_dir, paste0(study_name, "_eset_list.RData"))
annotations_path <- here::here(data_dir, "annotations", "full_chip_annotations.rds")
matrices_path <- here::here(data_dir, "expression_matrices", "global_cancer_matrix_maps.RData")

# ---- Ensure directories exist ----
dir.create(analysis_logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(aggregate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summaries_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(notebook_executed_dir, recursive = TRUE, showWarnings = FALSE)
