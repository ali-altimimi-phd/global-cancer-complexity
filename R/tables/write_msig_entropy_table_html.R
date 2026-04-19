write_msig_entropy_table_html <- function(comparison, entropy_df, output_dir = "quarto/resources/tables/msig") {
  fs::dir_create(output_dir)
  clean_name <- gsub("[^a-zA-Z0-9]", "_", tolower(comparison))
  out_file <- file.path(output_dir, paste0(clean_name, "_msig_entropy.html"))
  
  msig_entropy <- entropy_df %>%
    dplyr::filter(
      comparison == !!comparison,
      mode == "MSIG",
      p_perm <= 0.05
    ) %>%
    dplyr::mutate(
      group = dplyr::case_when(
        grepl("strongly anti-chaotic", spectral_direction, ignore.case = TRUE) ~ "Strongly Anti-Chaotic",
        grepl("mildly anti-chaotic",   spectral_direction, ignore.case = TRUE) ~ "Mildly Anti-Chaotic",
        grepl("neutral",               spectral_direction, ignore.case = TRUE) ~ "Neutral",
        grepl("mildly chaotic",        spectral_direction, ignore.case = TRUE) ~ "Mildly Chaotic",
        grepl("strongly chaotic",      spectral_direction, ignore.case = TRUE) ~ "Strongly Chaotic",
        TRUE ~ "Uncategorized"
      ),
      group = factor(
        group,
        levels = c(
          "Strongly Anti-Chaotic", "Mildly Anti-Chaotic", "Neutral",
          "Mildly Chaotic", "Strongly Chaotic", "Uncategorized"
        )
      )
    ) %>%
    dplyr::arrange(group, gene_set_name, p_perm)
  
  if (nrow(msig_entropy) == 0) {
    cat("No significant MSigDB gene sets found (spectral entropy).\n", file = out_file)
    return(invisible(NULL))
  }
  
  msig_entropy_split <- dplyr::group_split(msig_entropy, group)
  
  html_blocks <- lapply(msig_entropy_split, function(subgroup) {
    subgroup_title <- unique(subgroup$group)
    
    tbl <- subgroup %>%
      dplyr::select(GeneSet = gene_set_name, `p-value` = p_perm) %>%
      dplyr::arrange(GeneSet) %>%
      dplyr::mutate(`p-value` = sprintf("%.3f", `p-value`))
    
    table_html <- tbl %>%
      kableExtra::kable("html", align = c("l", "r")) %>%
      kableExtra::kable_styling(
        bootstrap_options = c("striped", "hover", "condensed"),
        full_width = TRUE
      )
    
    paste0(
      "<h4>", subgroup_title, "</h4>\n",
      as.character(table_html)
    )
  })
  
  full_html <- paste(html_blocks, collapse = "\n\n")
  writeLines(full_html, out_file)
  invisible(out_file)
}
