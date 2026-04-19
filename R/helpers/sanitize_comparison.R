# ------------------------------------------------------------------------------
# File: sanitize_comparison.R
# Purpose: Normalize comparison identifiers into filename-safe strings
# Role: Shared helper utility
# Pipeline: Shared
# Project: Cancer Complexity Analysis
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------

#' Sanitize comparison string for file and variable names
#'
#' Converts a comparison string (e.g., "BR/BRAD") to a lowercase,
#' filename-safe representation by replacing all non-alphanumeric
#' characters with underscores.
#'
#' @param x Character string to sanitize
#'
#' @return A sanitized, lowercase string safe for filenames or variable names
sanitize_comparison <- function(x) {
  gsub("[^A-Za-z0-9]", "_", tolower(x))
}