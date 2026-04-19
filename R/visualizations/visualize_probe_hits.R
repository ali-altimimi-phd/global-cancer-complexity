# R/visualizations/visualize_probe_hits.R

library(dplyr)
library(ggplot2)
library(glue)
library(forcats)

#' Plot filtered probe counts for a single comparison
#'
#' This function generates a bar plot showing the number of filtered probes
#' for each chip used in a specific comparison.
#'
#' @param df A data frame containing probe counts for a single comparison.
#'   Must include columns: `comparison`, `chip_label`, `filtered_probes`.
#'
#' @return A `ggplot` object (not rendered or saved; intended for export).
#' @export
plot_filtered_probe_counts <- function(df) {
  ggplot(df, aes(x = chip_label, y = filtered_probes, fill = chip_label)) +
    geom_col(width = 0.6, color = "black") +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = paste("Comparison:", unique(df$comparison)),
      x = "Chip",
      y = "Number of Filtered Probes"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")
}

#' Generate and save probe count plots for each comparison
#'
#' Loops over each unique comparison in the input data and saves a bar plot
#' of filtered probe counts across chips to PNG files in the output directory.
#'
#' @param probe_summary_df A data frame with columns: `comparison`, `chip_label`, and `filtered_probes`.
#' @param output_dir Directory where PNG files will be saved (default: `"quarto/resources/plots/probes"`).
#' @param plot_utils_path Path to the helper script containing `save_png_light()` (default: `"R/helpers/plot_utils.R"`).
#'
#' @return NULL. PNG files are saved to disk.
#' @export
plot_probe_counts_per_comparison <- function(probe_summary_df,
                                             output_dir = "quarto/resources/plots/probes",
                                             plot_utils_path = "R/helpers/plot_utils.R") {
  if (!file.exists(plot_utils_path)) {
    stop("plot_utils.R not found at: ", plot_utils_path)
  }
  source(plot_utils_path)  # for save_png_light()
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  comparisons <- unique(probe_summary_df$comparison)
  
  for (cmp in comparisons) {
    cmp_df <- dplyr::filter(probe_summary_df, comparison == cmp)
    if (nrow(cmp_df) == 0) next
    
    p <- plot_filtered_probe_counts(cmp_df)
    
    fname <- file.path(output_dir, paste0(sanitize_comparison(cmp), "_probe_counts.png"))
    save_png_light(p, filename = fname)
  }
}

#' Plot filtered probe counts grouped by cancer type
#'
#' Generates a grouped bar plot showing the number of filtered probes
#' per comparison across chips, grouped by cancer category.
#'
#' @param df A data frame with columns: `group`, `comparison`, `chip_label`, `filtered_probes`.
#'
#' @return A ggplot object (grouped bar plot).
#' @export
plot_probe_counts_by_group <- function(df) {
  df <- df %>%
    dplyr::mutate(
      comparison = forcats::fct_inorder(comparison),
      chip_label = factor(chip_label)
    )
  
  ggplot(df, aes(x = comparison, y = filtered_probes, fill = chip_label)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black") +
    scale_fill_brewer(palette = "Set2", name = "Chip") +
    facet_wrap(~ group, scales = "free_x", ncol = 1) +
    labs(
      title = "Filtered Probe Counts by Comparison and Chip",
      x = "Comparison",
      y = "Number of Filtered Probes"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      legend.position = "top"
    )
}

#' Save grouped probe count bar plots per cancer group
#'
#' Saves one PNG bar plot per group (e.g., carcinoma, leukemia), showing
#' side-by-side chip counts for each comparison.
#'
#' @param probe_df Data frame with columns: group, comparison, chip_label, filtered_probes
#' @param output_dir Directory to save plots (default: Quarto plot directory)
#' @export
save_probe_count_group_plots <- function(probe_df,
                                         output_dir = "quarto/resources/plots/probes") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  unique_groups <- unique(probe_df$group)
  
  for (grp in unique_groups) {
    grp_df <- dplyr::filter(probe_df, group == grp)
    if (nrow(grp_df) == 0) next
    
    p <- plot_probe_counts_by_group(grp_df)
    
    sanitized_group <- gsub("[^A-Za-z0-9_]", "_", tolower(grp))  # safe filename
    out_file <- file.path(output_dir, paste0(sanitized_group, "_grouped_probe_counts.png"))
    
    save_png_light(p, filename = out_file, width = 10, height = 6)
  }
}

#' Plot an UpSet plot showing overlap of filtered probes across chips
#'
#' @param chip_results A named list of chip result objects (e.g., list(hu35ksuba = res_hu35ksuba, hu6800 = res_hu6800))
#' @param chip_label_map Optional named character vector to relabel chips (e.g., c(hu35ksuba = "GPL98"))
#' @param output_file Full path to where the PNG plot should be saved
#' @param plot_utils_path Path to helper file containing save_png_light()
#'
#' @return NULL. Saves the UpSet plot to disk.
#' @export
plot_upset_probe_overlap <- function(chip_results,
                                     chip_label_map = NULL,
                                     output_file = "quarto/resources/plots/probes/upset_probe_overlap.png",
                                     plot_utils_path = "R/helpers/plot_utils.R") {
  stopifnot(requireNamespace("ComplexUpset", quietly = TRUE))
  
  if (!file.exists(plot_utils_path)) {
    stop("plot_utils.R not found at: ", plot_utils_path)
  }
  source(plot_utils_path)
  
  # Build union of all filtered probes
  all_probes <- Reduce(union, lapply(chip_results, function(res) res$filtered_probes))
  chip_names <- names(chip_results)
  
  # Build presence/absence matrix
  probe_matrix <- tibble::tibble(probe_id = all_probes)
  
  for (chip in chip_names) {
    chip_label <- if (!is.null(chip_label_map[[chip]])) chip_label_map[[chip]] else chip
    chip_set <- chip_results[[chip]]$filtered_probes
    probe_matrix[[chip_label]] <- all_probes %in% chip_set
  }
  
  # Generate UpSet plot
  p_upset <- ComplexUpset::upset(
    probe_matrix,
    intersect = colnames(probe_matrix)[-1],
    name = "Probe Presence",
    width_ratio = 0.2
  )
  
  # Save using standard light background
  save_png_light(p_upset, filename = output_file)
}
