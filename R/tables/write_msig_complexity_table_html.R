write_msig_complexity_table_html <- function(comparison, complexity_df, output_dir = "quarto/resources/tables/msig") {
  fs::dir_create(output_dir)
  clean_name <- gsub("[^a-zA-Z0-9]", "_", tolower(comparison))
  out_file <- file.path(output_dir, paste0(clean_name, "_msig_complexity.html"))
  
  msig_complex <- complexity_df %>%
    dplyr::filter(comparison == !!comparison, mode == "MSIG", p_perm <= 0.05) %>%
    dplyr::transmute(term = gene_set_name, p = p_perm, direction) %>%
    dplyr::distinct(term, .keep_all = TRUE) %>%
    dplyr::mutate(direction = factor(direction, levels = c("gained", "lost"))) %>%
    dplyr::arrange(direction, term)
  
  if (nrow(msig_complex) == 0) {
    writeLines("<p><em>No significant MSigDB gene sets found.</em></p>", out_file)
    return(invisible(out_file))
  }
  
  msig_complex_split <- dplyr::group_split(msig_complex, direction)
  
  html_blocks <- purrr::map_chr(msig_complex_split, function(subgroup) {
    subgroup_title <- stringr::str_to_title(unique(subgroup$direction))
    
    tbl <- subgroup %>%
      dplyr::select(GeneSet = term, `p-value` = p) %>%
      dplyr::arrange(GeneSet) %>%
      dplyr::mutate(`p-value` = sprintf("%.3f", `p-value`))
    
    html_table <- tbl %>%
      kableExtra::kable(format = "html", align = c("l", "r")) %>%
      kableExtra::kable_styling(
        bootstrap_options = c("striped", "hover", "condensed"),
        full_width = TRUE
      ) %>%
      as.character()
    
    paste0("<h4>", subgroup_title, " Complexity</h4>\n", html_table)
  })
  
  writeLines(html_blocks, out_file)
  invisible(out_file)
}
