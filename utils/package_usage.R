# ------------------------------------------------------------------------------
# File: package_usage.R
# Purpose: Utilities for scanning project files and summarizing package usage
# Role: Infrastructure utility
# Pipeline: Outside of pipelines
# Project: Cancer Complexity Analysis
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------

#' Package Usage Utilities
#'
#' Functions for scanning an R project and summarizing package usage across
#' R scripts, R Markdown, and Quarto documents.
#'
#' Supported patterns:
#' - library(pkg)
#' - require(pkg)
#' - pkg::fun
#'
#' Notes:
#' - This file defines reusable functions only; it does not execute work
#'   when sourced.
#' - Dynamic package loading (e.g., library(x) where x is a variable) is
#'   not detected.
#'

# ---- internal validation helpers --------------------------------------------

#' Validate required packages for dependency scanning
#'
#' @return invisible(TRUE)
#' @keywords internal
validate_package_usage_dependencies <- function() {
  required_pkgs <- c("stringr", "dplyr", "purrr", "tibble", "readr")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing_pkgs) > 0) {
    stop(
      paste0(
        "Missing required packages: ",
        paste(missing_pkgs, collapse = ", "),
        ". Please install them before running package usage scans."
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

#' Normalize file paths for regex matching
#'
#' @param x Character vector of file paths.
#' @return Character vector with forward slashes.
#' @keywords internal
normalize_scan_paths <- function(x) {
  gsub("\\\\", "/", x)
}

#' Determine whether a file should be scanned
#'
#' @param file Character scalar. File path.
#' @param exclude_patterns Character vector of regex patterns.
#' @return Logical scalar.
#' @keywords internal
should_scan_file <- function(file, exclude_patterns = NULL) {
  file_norm <- normalize_scan_paths(file)
  
  if (is.null(exclude_patterns) || length(exclude_patterns) == 0) {
    return(TRUE)
  }
  
  !any(stringr::str_detect(file_norm, exclude_patterns))
}

# ---- line-level extraction helpers ------------------------------------------

#' Extract packages loaded via library() or require()
#'
#' @param lines Character vector of lines from a file.
#' @param file Character scalar. Source file name.
#' @return Tibble with file, pkg, source, and matched text.
extract_library_require_calls <- function(lines, file = NA_character_) {
  tibble::tibble(
    file = file,
    text = lines
  ) %>%
    dplyr::mutate(
      pkg = stringr::str_match(
        text,
        "\\b(?:library|require)\\s*\\(\\s*['\"]?([A-Za-z][A-Za-z0-9._]*)['\"]?"
      )[, 2],
      source = "library_require"
    ) %>%
    dplyr::filter(!is.na(pkg)) %>%
    dplyr::select(file, pkg, source, text)
}

#' Extract packages referenced via pkg::fun syntax
#'
#' @param lines Character vector of lines from a file.
#' @param file Character scalar. Source file name.
#' @return Tibble with file, pkg, source, and matched text.
extract_namespace_calls <- function(lines, file = NA_character_) {
  matches <- stringr::str_extract_all(
    lines,
    "\\b([A-Za-z][A-Za-z0-9._]*)::"
  )
  
  pkg_vec <- unlist(matches, use.names = FALSE)
  pkg_vec <- gsub("::$", "", pkg_vec)
  
  if (length(pkg_vec) == 0) {
    return(
      tibble::tibble(
        file = character(),
        pkg = character(),
        source = character(),
        text = character()
      )
    )
  }
  
  tibble::tibble(
    file = file,
    pkg = pkg_vec,
    source = "namespace",
    text = NA_character_
  )
}

#' Extract package usage from one file's lines
#'
#' @param lines Character vector of lines.
#' @param file Character scalar. File path.
#' @return Tibble of detected package usage.
extract_package_usage_from_lines <- function(lines, file = NA_character_) {
  dplyr::bind_rows(
    extract_library_require_calls(lines, file = file),
    extract_namespace_calls(lines, file = file)
  )
}

# ---- file discovery ----------------------------------------------------------

#' List project files eligible for package scanning
#'
#' @param project_dir Character scalar. Root directory to scan.
#' @param include_extensions Character vector of file extensions without dots.
#' @param exclude_patterns Character vector of regex patterns to exclude.
#' @return Character vector of file paths.
list_package_scan_files <- function(
    project_dir = ".",
    include_extensions = c("R", "r", "Rmd", "rmd", "qmd"),
    exclude_patterns = c(
      "(^|/)renv/",
      "(^|/)docs/",
      "(^|/)_site/",
      "(^|/)site_libs/",
      "(^|/)\\.git/",
      "(^|/)\\.quarto/",
      "(^|/)archive/",
      "(^|/)old/"
    )
) {
  if (!dir.exists(project_dir)) {
    stop("project_dir does not exist: ", project_dir, call. = FALSE)
  }
  
  ext_pattern <- paste(include_extensions, collapse = "|")
  file_pattern <- paste0("\\.(", ext_pattern, ")$")
  
  files <- list.files(
    path = project_dir,
    pattern = file_pattern,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  files <- normalize_scan_paths(files)
  files[vapply(files, should_scan_file, logical(1), exclude_patterns = exclude_patterns)]
}

# ---- project scan ------------------------------------------------------------

#' Scan a project for package usage
#'
#' @param project_dir Character scalar. Root directory.
#' @param include_extensions Character vector of file extensions to scan.
#' @param exclude_patterns Character vector of regex patterns for files/folders
#'   to exclude.
#' @return A list with:
#' \describe{
#'   \item{raw}{Raw package usage records}
#'   \item{counts}{Counts by package}
#'   \item{counts_by_source}{Counts by package and source type}
#'   \item{file_presence}{Distinct package-file combinations}
#'   \item{files_scanned}{Vector of scanned files}
#' }
scan_project_packages <- function(
    project_dir = ".",
    include_extensions = c("R", "r", "Rmd", "rmd", "qmd"),
    exclude_patterns = c(
      "(^|/)renv/",
      "(^|/)docs/",
      "(^|/)_site/",
      "(^|/)site_libs/",
      "(^|/)\\.git/",
      "(^|/)\\.quarto/",
      "(^|/)archive/",
      "(^|/)old/"
    )
) {
  validate_package_usage_dependencies()
  
  files <- list_package_scan_files(
    project_dir = project_dir,
    include_extensions = include_extensions,
    exclude_patterns = exclude_patterns
  )
  
  if (length(files) == 0) {
    empty_raw <- tibble::tibble(
      file = character(),
      pkg = character(),
      source = character(),
      text = character()
    )
    
    return(list(
      raw = empty_raw,
      counts = tibble::tibble(pkg = character(), count = integer()),
      counts_by_source = tibble::tibble(
        pkg = character(),
        source = character(),
        count = integer()
      ),
      file_presence = tibble::tibble(
        pkg = character(),
        file = character(),
        source = character()
      ),
      files_scanned = character()
    ))
  }
  
  raw_results <- purrr::map_dfr(files, function(f) {
    lines <- readr::read_lines(f, progress = FALSE)
    extract_package_usage_from_lines(lines, file = f)
  })
  
  counts <- raw_results %>%
    dplyr::count(pkg, sort = TRUE, name = "count")
  
  counts_by_source <- raw_results %>%
    dplyr::count(pkg, source, sort = TRUE, name = "count")
  
  file_presence <- raw_results %>%
    dplyr::distinct(pkg, file, source)
  
  list(
    raw = raw_results,
    counts = counts,
    counts_by_source = counts_by_source,
    file_presence = file_presence,
    files_scanned = files
  )
}

# ---- convenience helpers -----------------------------------------------------

#' Identify candidate core packages for README documentation
#'
#' @param package_counts Tibble produced by scan_project_packages()$counts.
#' @param min_count Integer threshold for inclusion.
#' @return Tibble filtered to likely core dependencies.
identify_core_package_candidates <- function(package_counts, min_count = 5L) {
  if (!all(c("pkg", "count") %in% names(package_counts))) {
    stop("package_counts must contain columns 'pkg' and 'count'.", call. = FALSE)
  }
  
  package_counts %>%
    dplyr::filter(count >= min_count) %>%
    dplyr::arrange(dplyr::desc(count), pkg)
}

#' Write package usage outputs to disk
#'
#' @param scan_results List produced by scan_project_packages().
#' @param output_dir Character scalar. Directory for CSV outputs.
#' @param prefix Character scalar. Optional prefix for output filenames.
#' @return Invisible named character vector of written file paths.
write_package_usage_outputs <- function(
    scan_results,
    output_dir = ".",
    prefix = ""
) {
  required_names <- c("raw", "counts", "counts_by_source", "file_presence")
  missing_names <- setdiff(required_names, names(scan_results))
  
  if (length(missing_names) > 0) {
    stop(
      "scan_results is missing required elements: ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  prefix_str <- if (nzchar(prefix)) paste0(prefix, "_") else ""
  
  out_paths <- c(
    raw = file.path(output_dir, paste0(prefix_str, "package_usage_raw.csv")),
    counts = file.path(output_dir, paste0(prefix_str, "package_counts.csv")),
    counts_by_source = file.path(output_dir, paste0(prefix_str, "package_counts_by_source.csv")),
    file_presence = file.path(output_dir, paste0(prefix_str, "package_file_presence.csv"))
  )
  
  readr::write_csv(scan_results$raw, out_paths["raw"])
  readr::write_csv(scan_results$counts, out_paths["counts"])
  readr::write_csv(scan_results$counts_by_source, out_paths["counts_by_source"])
  readr::write_csv(scan_results$file_presence, out_paths["file_presence"])
  
  invisible(out_paths)
}