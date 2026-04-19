# entropy_plots.R
# -------------------------------------------------------------
# Visualization functions for Shannon entropy analysis
# -------------------------------------------------------------

source("R/plot_utils.R")
library(tibble)
library(ggplot2)

# 1. Bar plot for Normal vs. Cancer entropy with delta
plot_entropy_comparison <- function(tissue, entropy_n, entropy_c) {
  delta <- round(entropy_c - entropy_n, 2)
  
  make_ci_barplot(
    values = c(entropy_n, entropy_c),
    ci_lows = c(NA, NA),
    ci_highs = c(NA, NA),
    group_labels = c("Normal", "Cancer"),
    title = glue::glue("Shannon Entropy: {tissue} (Δ = {delta})"),
    ylab = "Entropy (bits)"
  )
}

# 2. Histogram of entropy values (optional, if you have per-sample entropy later)
plot_entropy_distribution <- function(entropies, group_labels, tissue) {
  df <- tibble(
    Entropy = entropies,
    Group = factor(group_labels, levels = unique(group_labels))
  )
  
  ggplot(df, aes(x = Entropy, fill = Group)) +
    geom_histogram(position = "identity", alpha = 0.6, bins = 30) +
    facet_wrap(~ Group, scales = "free_y") +
    theme_minimal() +
    ggtitle(paste("Entropy Distribution –", tissue)) +
    xlab("Shannon Entropy")
}
