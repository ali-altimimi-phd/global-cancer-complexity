# R/visualizations/generate_comparison_plots.R

#' Generate comparison-specific plots for a given comparison
#'
#' @param comparison Character name of the comparison (e.g., "PB/T-ALL")
#' @param complexity_df Filtered complexity results
#' @param entropy_df Filtered entropy results
#' @param plot_utils_path Path to plot_utils.R
#' @param output_dir Base directory for plots
#' @export
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
    dplyr::filter(comparison == !!comparison, mode == "GO", p_perm <= 0.05)
  
  if (nrow(go_comp) > 0) {
    p1 <- ggplot(go_comp, aes(x = -log10(p_perm), y = reorder(gene_set_name, p_perm), fill = direction)) +
      geom_col() +
      labs(title = paste("GO Complexity:", comparison), x = "-log10(p)", y = NULL) +
      theme_minimal()
    save_png_light(p1, file.path(plot_dirs$go, paste0(clean_id, "_go_complexity.png")))
  }
  
  ## ---- GO: Entropy ----
  go_ent <- entropy_df %>%
    dplyr::filter(comparison == !!comparison, mode == "GO", p_perm <= 0.05)
  
  if (nrow(go_ent) > 0) {
    p2 <- ggplot(go_ent, aes(x = -log10(p_perm), y = reorder(gene_set_name, p_perm), fill = spectral_direction)) +
      geom_col() +
      labs(title = paste("GO Entropy:", comparison), x = "-log10(p)", y = NULL) +
      theme_minimal()
    save_png_light(p2, file.path(plot_dirs$go, paste0(clean_id, "_go_entropy.png")))
  }
  
  ## ---- KEGG: Complexity ----
  kegg_comp <- complexity_df %>%
    dplyr::filter(comparison == !!comparison, mode == "KEGG", p_perm <= 0.05)
  
  if (nrow(kegg_comp) > 0) {
    p3 <- ggplot(kegg_comp, aes(x = -log10(p_perm), y = reorder(gene_set_name, p_perm), fill = direction)) +
      geom_col() +
      labs(title = paste("KEGG Complexity:", comparison), x = "-log10(p)", y = NULL) +
      theme_minimal()
    save_png_light(p3, file.path(plot_dirs$kegg, paste0(clean_id, "_kegg_complexity.png")))
  }
  
  ## ---- KEGG: Entropy ----
  kegg_ent <- entropy_df %>%
    dplyr::filter(comparison == !!comparison, mode == "KEGG", p_perm <= 0.05)
  
  if (nrow(kegg_ent) > 0) {
    p4 <- ggplot(kegg_ent, aes(x = -log10(p_perm), y = reorder(gene_set_name, p_perm), fill = spectral_direction)) +
      geom_col() +
      labs(title = paste("KEGG Entropy:", comparison), x = "-log10(p)", y = NULL) +
      theme_minimal()
    save_png_light(p4, file.path(plot_dirs$kegg, paste0(clean_id, "_kegg_entropy.png")))
  }
  
  ## ---- KEGG Cancer Pathways ----
  kegg_cancer_terms <- bind_rows(
    complexity_df %>% filter(comparison == !!comparison, mode == "KEGG", kegg_cancer == TRUE),
    entropy_df %>% filter(comparison == !!comparison, mode == "KEGG", kegg_cancer == TRUE)
  ) %>%
    dplyr::select(gene_set_name, p_perm, direction = spectral_direction, source = chip) %>%
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
}
