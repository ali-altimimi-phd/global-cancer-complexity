# R/config/global_cancer/analysis_config.R

# ---- Study Identifier ----
# Used to define directory structure for multiple data sources
study_name <- "global_cancer"  # used to define directory structure

# ---- Analysis Pipeline Toggles ----
build_matrix_maps <- FALSE # Toggle to enable or disable pipeline stage
run_filtering <- TRUE
run_pairwise <- TRUE
run_aggregator <- TRUE
run_comparison_summary <- TRUE

# ---- Filtering Parameters ----
filter_method <- "limma"  # "limma" or "variance"
logfc_cutoff <- 0.33
pval_cutoff <- 0.05
var_threshold <- 0.75
top_n <- 3000

# ---- Pairwise Comparison Settings ----
# Chip platforms to analyze
chips <- c("hu35ksuba", "hu6800")

geo_chip_map <- list(hu35ksuba = "GPL98", hu6800 = "GPL80")

# Engines to run: supports 'complexity' and/or 'entropy'
engines <- c("complexity", "entropy")

# ---- Gene Set Filter ----
gene_set_mode <- "ALL"           # Options: "ALL", "GO_BP", "GO_MF", "KEGG", "MSIGDB"
# gene_set_mode <- "GO_BP"
# gene_set_mode <- "GO_MF"
# gene_set_mode <- "KEGG"          
# gene_set_mode <- "MSIGDB"
quantile_cutoff <- 0.75          # Dynamic threshold (used with GO/KEGG/MSIGDB)
min_probes <- 5                  # Absolute minimum probes per gene set

# ---- Local Paths ----
logs_dir      <- here::here("output", study_name, "logs")
data_dir      <- here::here("output", study_name, "RData")
filtered_probes_dir <- here::here(data_dir, "filtered_probes")
aggregate_dir <- here::here(data_dir, "aggregated")
summaries_dir  <- here::here(data_dir, "summaries")
reports_dir    <- here::here("output", study_name, "reports")

# ---- Output File Paths ----
analysis_pipeline_logfile <- here::here(logs_dir, "analysis", "analysis_pipeline_log.txt")
eset_path        <- here::here(data_dir, paste0(study_name, "_eset_list.RData"))
annotations_path <- here::here(data_dir, "annotations", "full_chip_annotations.rds")
matrices_path <- here::here(data_dir, "expression_matrices", "global_cancer_matrix_maps.RData")

# ---- Ensure directories exist ----
dir.create(dirname(analysis_pipeline_logfile), recursive = TRUE, showWarnings = FALSE)
dir.create(aggregate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summaries_dir, recursive = TRUE, showWarnings = FALSE)
