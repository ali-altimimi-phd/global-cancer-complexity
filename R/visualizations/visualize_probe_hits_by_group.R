# R/visualizations/visualize_probe_hits_by_group.R
# --------------------------------------------------

library(dplyr)
library(ggplot2)
library(forcats)
library(glue)
library(ComplexUpset)

#' Plot filtered probe counts by group
#'
#' @param df A data frame with: group, comparison, chip_label, filtered_probes
#'
#' @return A ggplot object
#' @export
plot_probe_counts_by_group <- function(df) {
  df <- df %>%
    mutate(
      comparison = forcats::fct_inorder(comparison),
      chip_label = factor(chip_label)
    )
  
  ggplot(df, aes(x = comparison, y = filtered_probes, fill = chip_label)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black") +
    scale_fill_brewer(palette = "Set2", name = "Chip") +
    labs(
      title = paste("Filtered Probe Counts -", unique(df$group)),
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

#' Plot UpSet plot of filtered probe overlaps (by comparison)
#'
#' @param df A data frame with columns: group, comparison, filtered_probes (list-column)
#' @param min_frac Minimum fraction of comparisons a probe must appear in (default: 0.3 if n >= 4)
#'
#' @return A ComplexUpset ggplot object
#' @export
plot_upset_probe_overlap_by_group <- function(df, min_frac = 0.3) {
  stopifnot("filtered_probes" %in% colnames(df))
  
  comparisons <- df$comparison
  all_probes <- unique(unlist(df$filtered_probes))
  
  # Build logical presence/absence matrix
  upset_matrix <- tibble::tibble(probe_id = all_probes)
  for (cmp in comparisons) {
    probes <- df$filtered_probes[df$comparison == cmp][[1]]
    upset_matrix[[cmp]] <- all_probes %in% probes
  }
  
  if (length(comparisons) < 4) {
    return(plot_simple_upset(upset_matrix))
  } else {
    # Try advanced plot, fallback to simple if it fails
    tryCatch(
      {
        plot_filtered_upset(upset_matrix, comparisons, min_frac)
      },
      error = function(e) {
        message("⚠️ Falling back to simple UpSet plot due to error: ", e$message)
        plot_simple_upset(upset_matrix)
      }
    )
  }
}

# ---- Subfunction A: simple plot ----
plot_simple_upset <- function(upset_matrix) {
  ComplexUpset::upset(
    upset_matrix,
    intersect = colnames(upset_matrix)[-1],
    name = "Comparison",
    width_ratio = 0.2
  )
}

# ---- Subfunction B: filtered with custom label ----
plot_filtered_upset <- function(upset_matrix, comparisons, min_frac) {
  row_counts <- rowSums(upset_matrix[, -1])
  min_required <- max(3, floor(length(comparisons) * min_frac))
  upset_matrix <- upset_matrix[row_counts >= min_required, , drop = FALSE]
  
  label <- glue::glue("intersection size (≥ {min_required} comps)")
  
  ComplexUpset::upset(
    upset_matrix,
    intersect = colnames(upset_matrix)[-1],
    name = "Comparison",
    width_ratio = 0.2,
    base_annotations = setNames(
      list(
        ComplexUpset::intersection_size(
          counts = TRUE
        )
      ),
      label  # <- interpolated label as the name
    )
  )
}

#' Save grouped bar plots and UpSet plots for each cancer group
#'
#' @param probe_df A data frame with: group, comparison, chip_label, filtered_probes
#' @param output_dir Output directory to save PNGs (default: quarto probe group dir)
#' @param plot_utils_path Path to helper script (must include save_png_light())
#'
#' @return NULL. Plots are saved to disk.
#' @export
save_probe_group_visualizations <- function(probe_df,
                                            output_dir = "quarto/resources/plots/probes/groups",
                                            plot_utils_path = "R/helpers/plot_utils.R") {
  if (!file.exists(plot_utils_path)) {
    stop("plot_utils.R not found at: ", plot_utils_path)
  }
  source(plot_utils_path)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  unique_groups <- unique(probe_df$group)
  
  for (grp in unique_groups) {
    grp_df <- filter(probe_df, group == grp)
    if (nrow(grp_df) == 0) next
    
    # Save bar plot
    p_bar <- plot_probe_counts_by_group(grp_df)
    bar_file <- file.path(output_dir, paste0(gsub("[^A-Za-z0-9_]", "_", tolower(grp)), "_probe_counts.png"))
    save_png_light(p_bar, filename = bar_file, width = 10, height = 6)
    
    # Only generate UpSet if we have filtered probes for 2+ comparisons
    if (length(unique(grp_df$comparison)) >= 2) {
      try({
        # Set quantile threshold here, e.g., Q10
        p_upset <- plot_upset_probe_overlap_by_group(grp_df)
        upset_file <- file.path(output_dir, paste0(gsub("[^A-Za-z0-9_]", "_", tolower(grp)), "_probe_overlap_upset.png"))
        save_png_light(p_upset, filename = upset_file, width = 8, height = 5)
      }, silent = TRUE)
    }
  }
}
