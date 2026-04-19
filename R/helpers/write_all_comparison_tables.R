#' Write All HTML Tables for a Comparison
#'
#' Writes GO, KEGG, and MSigDB tables (complexity and entropy) for a given comparison.
#' 
#' Used in `run_all_reports()` to reduce clutter and increase modularity.
#'
#' @param comparison Comparison code (e.g., "BR/BRAD")
#' @param complexity_df Complexity results data frame
#' @param entropy_df Entropy results data frame
#' @param output_dir Output directory for HTML files
#' @export
write_all_comparison_tables <- function(comparison,
                                        complexity_df,
                                        entropy_df,
                                        output_dir) {

  source(here::here("R/tables/write_go_complexity_table_html.R"))
  source(here::here("R/tables/write_go_entropy_table_html.R"))
  source(here::here("R/tables/write_kegg_complexity_table_html.R"))
  source(here::here("R/tables/write_kegg_entropy_table_html.R"))
  source(here::here("R/tables/write_msig_complexity_table_html.R"))
  source(here::here("R/tables/write_msig_entropy_table_html.R"))

  write_go_complexity_table_html(
    comparison    = comparison,
    complexity_df = complexity_df,
    output_dir    = output_dir
  )
  
  write_go_entropy_table_html(comparison    = comparison,
                              entropy_df    = entropy_df,
                              output_dir    = output_dir)
  
  write_kegg_complexity_table_html(
    comparison    = comparison,
    complexity_df = complexity_df,
    output_dir    = output_dir
  )
  
  write_kegg_entropy_table_html(comparison    = comparison,
                                entropy_df    = entropy_df,
                                output_dir    = output_dir)
  
  write_msig_complexity_table_html(
    comparison    = comparison,
    complexity_df = complexity_df,
    output_dir    = output_dir
  )
  
  write_msig_entropy_table_html(comparison    = comparison,
                                entropy_df    = entropy_df,
                                output_dir    = output_dir)
}
