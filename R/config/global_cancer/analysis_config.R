# ------------------------------------------------------------------------------
# File: analysis_config.R
# Purpose: Configuration for the global cancer analysis pipeline
# Role: Defines toggles, filtering parameters, latent-space settings, paths,
#       output locations, and logging destinations.
# Pipeline: Analysis
# Project: Cancer Complexity Analysis
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------

# ==============================================================================
# 1. Study identifier
# ==============================================================================

study_name <- "global_cancer"


# ==============================================================================
# 2. Main analysis pipeline toggles
# ==============================================================================

build_matrix_maps      <- FALSE
run_filtering          <- FALSE
run_pairwise           <- FALSE
run_aggregator         <- FALSE
run_comparison_summary <- FALSE


# ==============================================================================
# 3. Probe-selection and filtering parameters
# ==============================================================================

# Primary probe-selection method:
#   "limma"     = differential-expression-based probe selection
#   "variance" = variability-based probe selection
filter_method <- "limma"

# ---- Limma mode ----
# Used only when filter_method == "limma"
logfc_cutoff <- 0.33
pval_cutoff  <- 0.05

# ---- Variance mode ----
# Used only when filter_method == "variance"
#
# variance_selection_mode:
#   "top_n"     = retain a fixed number of most variable probes
#   "threshold" = retain probes above a variability quantile cutoff
variance_selection_mode <- "top_n"

# Used only when variance_selection_mode == "top_n"
analysis_top_n <- 3000

# Used only when variance_selection_mode == "threshold"
var_threshold <- 0.75


# ==============================================================================
# 4. Pairwise comparison and gene-set settings
# ==============================================================================

chips <- c("hu35ksuba", "hu6800")

geo_chip_map <- list(
  hu35ksuba = "GPL98",
  hu6800    = "GPL80"
)

engines <- c("complexity", "entropy")

# Gene-set modes:
#   "FULL"   = probe-level / full selected feature space
#   "GO_BP"  = Gene Ontology biological process
#   "GO_MF"  = Gene Ontology molecular function
#   "KEGG"   = KEGG pathway gene sets
#   "MSIGDB" = MSigDB Hallmark gene sets
# gene_set_modes <- c("FULL", "GO_BP", "GO_MF", "KEGG", "MSIGDB")
gene_set_modes <- c("GO_MF", "KEGG", "MSIGDB")

# Gene-set-level filtering parameters
quantile_cutoff <- 0.75
min_probes      <- 5


# ==============================================================================
# 5. Core output and logging paths
# ==============================================================================

logs_dir          <- here::here("output", study_name, "logs")
analysis_logs_dir <- file.path(logs_dir, "analysis")

data_dir            <- here::here("output", study_name, "RData")
filtered_probes_dir <- file.path(data_dir, "filtered_probes")
aggregate_dir       <- file.path(data_dir, "aggregated")
summaries_dir       <- file.path(data_dir, "summaries")
reports_dir         <- here::here("output", study_name, "reports")

analysis_pipeline_logfile <- file.path(
  analysis_logs_dir,
  "analysis_pipeline_log.txt"
)

eset_path <- file.path(
  data_dir,
  paste0(study_name, "_eset_list.RData")
)

annotations_path <- file.path(
  data_dir,
  "annotations",
  "full_chip_annotations.rds"
)

matrices_path <- file.path(
  data_dir,
  "expression_matrices",
  "global_cancer_matrix_maps.RData"
)


# ==============================================================================
# 6. Latent-space project paths and Python environment
# ==============================================================================

# Project-level latent-space source directory
latent_project_dir <- here::here("projects", "cancer-latent-space")
latent_script_dir  <- file.path(latent_project_dir, "scripts")
latent_notebook_dir <- file.path(latent_project_dir, "notebooks")

# Python executable used by both script execution and Jupyter nbconvert.
# Can be overridden externally by setting LATENT_PYTHON_EXE.
latent_python_exe <- Sys.getenv(
  "LATENT_PYTHON_EXE",
  unset = "C:/Users/drziy/anaconda3/envs/ml/python.exe"
)


# ==============================================================================
# 7. Latent-space Python script execution
# ==============================================================================

# Runs production Python scripts replacing notebooks 00--03:
#   run_latent_preparation.py
#   run_latent_training.py
run_latent_python_scripts <- FALSE

latent_python_scripts <- c(
  "run_latent_preparation.py",
  "run_latent_training.py"
)

# Latent-space preprocessing/model parameters
latent_chip_id        <- "hu35ksuba"
latent_filter_method  <- "variance"
latent_top_n          <- 3000

latent_dim            <- 10
latent_epochs         <- 50
latent_batch_size     <- 32
latent_beta           <- 0.01
latent_learning_rate  <- 1e-3


# ==============================================================================
# 8. Latent-space notebook execution
# ==============================================================================

# Runs remaining analysis notebooks after Python scripts have produced latent data.
# Notebooks 00--03 are no longer listed here because they are replaced by scripts.
run_latent_notebooks <- FALSE

execute_notebooks_inplace <- FALSE
notebook_timeout_seconds  <- -1

jupyter_command     <- latent_python_exe
jupyter_prefix_args <- c("-m", "jupyter")

notebook_executed_dir <- here::here(
  "output",
  study_name,
  "notebooks",
  "executed"
)

latent_notebooks <- c(
  "04_latent_space_analysis_standardized.ipynb",
  "05_latent_complexity_analysis_standardized.ipynb",
  "06_between_cancer_geometry_in_latent_space_enhanced.ipynb"
)


# ==============================================================================
# 9. Latent-space postprocessing
# ==============================================================================

# Builds integrated latent + legacy outputs:
#   master_combined_metrics
#   latent_taxonomy_table
#   regime summaries
#   separation_vs_complexity plot
run_latent_postprocessing <- TRUE

# Output-level latent-space directory.
# Deliberately distinct from latent_project_dir, which points to source code.
latent_output_dir <- here::here("output", study_name)

latent_postprocessing_table_dir <- file.path(
  latent_output_dir,
  "tables",
  "latent"
)

latent_postprocessing_plot_dir <- file.path(
  latent_output_dir,
  "plots",
  "latent"
)

latent_postprocessing_chip_id       <- "hu35ksuba"
latent_postprocessing_mode          <- "ALL"
latent_postprocessing_gene_set_name <- "ALL"

# Complexity/structure delta column used for latent taxonomy classification.
# This should be explicit rather than auto-detected, because both legacy-style
# and report-style columns may exist after aggregation.
latent_postprocessing_complexity_delta_col <- "pr_delta"

# ==============================================================================
# 10. Ensure required directories exist
# ==============================================================================

dir.create(analysis_logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(aggregate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summaries_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(notebook_executed_dir, recursive = TRUE, showWarnings = FALSE)


dir.create(
  latent_postprocessing_table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  latent_postprocessing_plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)