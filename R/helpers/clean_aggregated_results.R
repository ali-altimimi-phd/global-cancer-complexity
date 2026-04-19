#' Clean and standardize entropy/complexity summary data
#'
#' @param df A tibble from either entropy or complexity analysis.
#' @param drop_engine Logical, whether to drop the 'engine' column (default = TRUE)
#'
#' @return Cleaned tibble
clean_aggregated_results <- function(df, drop_engine = TRUE) {
  # --- Drop or rename redundant columns ---
  cols_to_drop <- c("Chip", "gene_set")
  if (drop_engine) cols_to_drop <- c(cols_to_drop, "engine")
  
  df <- df |>
    dplyr::select(-dplyr::any_of(cols_to_drop)) |>
    dplyr::rename_with(tolower) |>
    dplyr::rename(number_of_probes = filter) |>
    dplyr::relocate(mode, .before = number_of_probes)
  
  # --- Drop redundant lookup columns if readable names are present ---
  if ("gene_set_name" %in% names(df)) {
    df <- df |>
      dplyr::select(-dplyr::any_of(c("go_name", "path_name")))
  }
  
  # --- Reorder key columns ---
  desired_order <- c(
    "comparison", "chip", "mode", "number_of_probes",
    "gene_set_normalized", "gene_set_name", "go_ontology", "p_perm"
  )
  existing_order <- intersect(desired_order, names(df))
  rest <- setdiff(names(df), existing_order)
  
  df <- df |>
    dplyr::select(dplyr::all_of(c(existing_order, rest)))
  
  return(df)
}