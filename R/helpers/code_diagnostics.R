# code_diagnostics.R
# -------------------------------------------------------------
# Utility functions for searching and auditing codebase files
# Used to find patterns such as "annotations" across scripts
# -------------------------------------------------------------

#' Search R code files for a pattern
#'
#' Recursively searches all `.R` files in the "R/" directory for lines
#' containing a specified string or regex pattern.
#'
#' @param pattern String or regex to search for (e.g., "annotations").
#' @param path Base folder to start from (default = "R/").
#' @param ignore_case Logical; if TRUE, performs case-insensitive search.
#' @return A character vector of file paths containing the pattern.
#' @examples
#' # From R console or script:
#' search_r_files_for_pattern("annotations")
#' search_r_files_for_pattern("gene_set_filter", ignore_case = TRUE)
#'
#' @export
search_r_files_for_pattern <- function(pattern,
                                       path = "R/",
                                       ignore_case = FALSE) {
  # List all R files recursively
  r_files <- list.files(path, recursive = TRUE, full.names = TRUE, pattern = "\\.R$")
  
  matching_files <- character()
  
  for (f in r_files) {
    lines <- tryCatch(readLines(f, warn = FALSE), error = function(e) NULL)
    if (!is.null(lines)) {
      match_found <- any(grepl(pattern, lines, ignore.case = ignore_case))
      if (match_found) matching_files <- c(matching_files, f)
    }
  }
  
  return(matching_files)
}

# -------------------------------------------------------------
# Example usage (interactive only, not run during pipeline):
# -------------------------------------------------------------
# To search for files referencing 'annotations':
# source("R/helpers/code_diagnostics.R")
# matching <- search_r_files_for_pattern("annotations")
# cat(matching, sep = "\n")

# Add the ignore_case = TRUE argument if you're unsure of capitalization:
# search_r_files_for_pattern("gene_set", ignore_case = TRUE)
