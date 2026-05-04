# ------------------------------------------------------------------------------
# File: generate_comparison_plots.R
# Purpose: Generate comparison-level diagnostic plots for GO, KEGG, KEGG cancer
#   pathways, and filtered-probe summaries for Quarto reporting.
# Role: Helper (comparison plot generator)
# Pipeline: Reporting
# Project: Global Cancer Complexity
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------

#' Generate Comparison-Level Reporting Plots
#'
#' Generates comparison-specific plots summarizing significant GO terms, KEGG
#' pathways, KEGG cancer pathways, and filtered-probe counts for inclusion in
#' Quarto-generated reports.
#'
#' @param comparison Character string identifying the comparison
#'   (e.g., "PB/T-ALL").
#' @param complexity_df Data frame containing complexity results.
#' @param entropy_df Data frame containing entropy results.
#' @param plot_utils_path Path to plotting utility functions.
#' @param output_dir Base directory where plot files will be written.
#'
#' @return Invisibly returns \code{NULL}; writes PNG plots to disk.

generate_comparison_plots <- function(comparison,
                                      complexity_df,
                                      entropy_df,
                                      plot_utils_path = "R/helpers/plot_utils.R",
                                      output_dir = "quarto/resources/plots") {
  source(plot_utils_path)
  
  clean_id <- gsub("[^a-zA-Z0-9]", "_", tolower(comparison))
  plot_dirs <- list(
    go = file.path(output_dir, "go"),
    kegg = file.path(output_dir, "kegg"),
    kegg_cancer = file.path(output_dir, "kegg_cancer"),
    probe = file.path(output_dir, "probes")
  )
  lapply(plot_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  
  ## ---- GO: Complexity ----
  go_comp <- complexity_df %>%
    dplyr::filter(comparison == !!comparison, mode == "GO", p_perm <= 0.05) %>%
    valid_plot_rows()  
    
  if (nrow(go_comp) > 0) {
    p1 <- ggplot(go_comp, aes(x = -log10(p_perm), y = reorder(gene_set_name, p_perm), fill = direction)) +
      geom_col() +
      labs(title = paste("GO Complexity:", comparison), x = "-log10(p)", y = NULL) +
      theme_minimal()
    save_png_light(p1, file.path(plot_dirs$go, paste0(clean_id, "_go_complexity.png")))
  }
  
  ## ---- GO: Entropy ----
  go_ent <- entropy_df %>%
    dplyr::filter(comparison == !!comparison, mode == "GO", p_perm <= 0.05) %>%
    valid_plot_rows()  
  
  if (nrow(go_ent) > 0) {
    p2 <- ggplot(go_ent, aes(x = -log10(p_perm), y = reorder(gene_set_name, p_perm), fill = spectral_direction)) +
      geom_col() +
      labs(title = paste("GO Entropy:", comparison), x = "-log10(p)", y = NULL) +
      theme_minimal()
    save_png_light(p2, file.path(plot_dirs$go, paste0(clean_id, "_go_entropy.png")))
  }
  
  ## ---- KEGG: Complexity ----
  kegg_comp <- complexity_df %>%
    dplyr::filter(comparison == !!comparison, mode == "KEGG", p_perm <= 0.05) %>%
    valid_plot_rows()
  
  if (nrow(kegg_comp) > 0) {
    p3 <- ggplot(kegg_comp, aes(x = -log10(p_perm), y = reorder(gene_set_name, p_perm), fill = direction)) +
      geom_col() +
      labs(title = paste("KEGG Complexity:", comparison), x = "-log10(p)", y = NULL) +
      theme_minimal()
    save_png_light(p3, file.path(plot_dirs$kegg, paste0(clean_id, "_kegg_complexity.png")))
  }
  
  ## ---- KEGG: Entropy ----
  kegg_ent <- entropy_df %>%
    dplyr::filter(comparison == !!comparison, mode == "KEGG", p_perm <= 0.05) %>%
    valid_plot_rows()
  
  if (nrow(kegg_ent) > 0) {
    p4 <- ggplot(kegg_ent, aes(x = -log10(p_perm), y = reorder(gene_set_name, p_perm), fill = spectral_direction)) +
      geom_col() +
      labs(title = paste("KEGG Entropy:", comparison), x = "-log10(p)", y = NULL) +
      theme_minimal()
    save_png_light(p4, file.path(plot_dirs$kegg, paste0(clean_id, "_kegg_entropy.png")))
  }
  
  ## ---- KEGG Cancer Pathways ----
  kegg_cancer_complexity <- complexity_df %>%
    dplyr::filter(
      comparison == !!comparison,
      mode == "KEGG",
      kegg_cancer == TRUE
    ) %>%
    dplyr::transmute(
      gene_set_name = gene_set_name,
      p_perm = p_perm,
      direction = direction,
      source = paste0(chip, " / complexity")
    )
  
  kegg_cancer_entropy <- entropy_df %>%
    dplyr::filter(
      comparison == !!comparison,
      mode == "KEGG",
      kegg_cancer == TRUE
    ) %>%
    dplyr::transmute(
      gene_set_name = gene_set_name,
      p_perm = p_perm,
      direction = spectral_direction,
      source = paste0(chip, " / entropy")
    )
  
  kegg_cancer_terms <- dplyr::bind_rows(
    kegg_cancer_complexity,
    kegg_cancer_entropy
  ) %>%
    valid_plot_rows() %>%
    dplyr::distinct()  
    
  if (nrow(kegg_cancer_terms) > 0) {
    p5 <- ggplot(kegg_cancer_terms, aes(x = -log10(p_perm), y = reorder(gene_set_name, p_perm), fill = direction)) +
      geom_col() +
      labs(title = paste("KEGG Cancer Pathways:", comparison), x = "-log10(p)", y = NULL) +
      theme_minimal()
    save_png_light(p5, file.path(plot_dirs$kegg_cancer, paste0(clean_id, "_kegg_cancer_terms.png")))
  }
  
  ## ---- Probe Count ----
  probe_df <- data.frame(
    Chip = c("hu35ksuba", "hu6800"),
    Count = c(
      sum(complexity_df$comparison == comparison & complexity_df$chip == "hu35ksuba"),
      sum(complexity_df$comparison == comparison & complexity_df$chip == "hu6800")
    )
  )
  
  p6 <- ggplot(probe_df, aes(x = Chip, y = Count)) +
    geom_col(fill = "steelblue") +
    theme_minimal() +
    labs(title = paste("Filtered Probes:", comparison), y = "Probe Count")
  
  save_png_light(p6, file.path(plot_dirs$probe, paste0(clean_id, "_probe_counts.png")))
  
  invisible(NULL)
}
