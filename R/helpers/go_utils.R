#' Count GO terms by ontology (BP, MF, CC) using internal annotations
#'
#' @param df A data frame containing a `go_ontology` column.
#' @param mode_filter Optional mode to restrict to (e.g. "GO").
#'
#' @return A tibble with counts of GO terms per ontology.
#' @export
get_internal_go_ontology_counts <- function(df, mode_filter = NULL) {
  if (!"go_ontology" %in% colnames(df)) {
    stop("Data frame must contain a 'go_ontology' column.")
  }
  
  df <- df |> dplyr::filter(!is.na(go_ontology))
  
  if (!is.null(mode_filter) && "mode" %in% names(df)) {
    df <- df |> dplyr::filter(mode == mode_filter)
  }
  
  df |> dplyr::count(go_ontology, name = "n_terms") |> dplyr::arrange(desc(n_terms))
}
