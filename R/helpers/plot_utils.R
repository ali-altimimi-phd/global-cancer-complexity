#' plot_utils.R
#' -------------------------------------------------------------
#' General-purpose plotting utilities used across analyses
#' -------------------------------------------------------------

#' @author Ali M. Al-Timimi

library(dplyr)
library(ggplot2)
library(tibble)


#' Filter Rows with Valid P-Values for Plotting
#'
#' Removes rows with missing, non-finite, or non-positive permutation p-values
#' before computing transformations such as \code{-log10(p_perm)}.
#'
#' @param df Data frame containing a \code{p_perm} column.
#'
#' @return Filtered data frame containing only rows with valid positive p-values.
#' @export
valid_plot_rows <- function(df) {
  df |>
    dplyr::filter(
      !is.na(p_perm),
      is.finite(p_perm),
      p_perm > 0
    )
}

#' Create a bar plot with confidence intervals
#'
#' @param values Numeric vector of means.
#' @param ci_lows Numeric vector of lower CI bounds.
#' @param ci_highs Numeric vector of upper CI bounds.
#' @param group_labels Character vector of group labels.
#' @param title Plot title.
#' @param ylab Y-axis label.
#' @return A ggplot object.
#' @export
make_ci_barplot <- function(values, ci_lows, ci_highs, group_labels, title, ylab) {
  tibble(
    Group = group_labels,
    Mean = values,
    CI_Low = ci_lows,
    CI_High = ci_highs
  ) |>
    ggplot(aes(x = Group, y = Mean)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.2) +
    ggtitle(title) +
    ylab(ylab) +
    theme_minimal()
}

#' Create a histogram showing permutation test results
#'
#' @param null_dist Numeric vector of permuted values.
#' @param observed_value Numeric observed value.
#' @param title Plot title.
#' @param xlab X-axis label.
#' @return A ggplot object.
#' @export
make_permutation_plot <- function(null_dist, observed_value, title, xlab) {
  tibble(Permuted = null_dist) |>
    ggplot(aes(x = Permuted)) +
    geom_histogram(bins = 40, fill = "gray", color = "black") +
    geom_vline(xintercept = observed_value, color = "red", linetype = "dashed", linewidth = 1) +
    ggtitle(title) +
    xlab(xlab) +
    theme_minimal()
}

#' Create a violin plot with jittered points
#'
#' @param values Numeric vector of values.
#' @param groups Factor or character vector of group labels.
#' @param title Plot title.
#' @param ylab Y-axis label.
#' @return A ggplot object.
#' @export
make_violin_jitter_plot <- function(values, groups, title, ylab) {
  df <- tibble(Value = values, Group = groups)
  
  ggplot(df, aes(x = Group, y = Value, fill = Group)) +
    geom_violin(trim = FALSE, color = "black") +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1.5) +
    ggtitle(title) +
    ylab(ylab) +
    theme_minimal()
}

#' Save a ggplot to PNG with light background
#'
#' @param plot ggplot object to save.
#' @param filename Output file name (with .png extension).
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Resolution in dots per inch.
#' @return NULL. File is written to disk.
#' @export
save_png_light <- function(plot, filename, width = 7, height = 5, dpi = 300) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}
