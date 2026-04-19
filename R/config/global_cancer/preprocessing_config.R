# R/config/global_cancer/preprocessing_config.R

# ---- Study Identifier ----
# Used to define directory structure for multiple data sources
study_name <- "global_cancer"

# ---- Preprocessing Pipeline Toggles ----
download_enabled      <- FALSE
build_esets           <- FALSE
process_metadata      <- FALSE
run_annotation        <- TRUE

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
# These are specific to global cancer data set
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
preprocess_pipeline_logfile <- here::here(logs_dir, "preprocess", "preprocess_pipeline_log.txt")
eset_path        <- here::here(data_dir, paste0(study_name, "_eset_list.RData"))
annotations_path <- here::here(data_dir, "annotations", "full_chip_annotations.rds")

# ---- Ensure Directories Exist ----
dir.create(dirname(preprocess_pipeline_logfile), recursive = TRUE, showWarnings = FALSE)
dir.create(local_cel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(local_geo_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(eset_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(annotations_path), recursive = TRUE, showWarnings = FALSE)
