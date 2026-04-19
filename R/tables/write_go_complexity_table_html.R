# File: R/tables/write_go_complexity_table_html.R

write_go_complexity_table_html <- function(comparison, complexity_df, output_dir = "quarto/resources/tables/go") {
  fs::dir_create(output_dir)
  clean_name <- gsub("[^a-zA-Z0-9]", "_", tolower(comparison))
  out_file <- file.path(output_dir, paste0(clean_name, "_go_complexity.html"))
  
  go_complex <- complexity_df %>%
    dplyr::filter(comparison == !!comparison, mode == "GO", p_perm <= 0.05) %>%
    dplyr::mutate(direction = factor(direction, levels = c("gained", "lost"))) %>%
    dplyr::arrange(direction, gene_set_name, p_perm)
  
  if (nrow(go_complex) == 0) {
    writeLines("<p><em>No significant GO terms found.</em></p>", out_file)
    return(invisible(out_file))
  }
  
  go_complex_split <- dplyr::group_split(go_complex, direction)
  
  html_blocks <- purrr::map_chr(go_complex_split, function(subgroup) {
    subgroup_title <- stringr::str_to_title(unique(subgroup$direction))
    
    tbl <- subgroup %>%
      dplyr::select(GeneSet = gene_set_name, `p-value` = p_perm) %>%
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
