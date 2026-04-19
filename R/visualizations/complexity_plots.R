# complexity_plots.R
# -------------------------------------------------------------
# Wrappers around shared plotting utils for complexity output
# -------------------------------------------------------------

library(here)

source(here("R/plot_utils.R"))

# 1. Plot CI barplot for SVD kappa -----------------------------

plot_kappa_ci <- function(tissue, res) {
  make_ci_barplot(
    values = c(res$`SVD kappa Normal`, res$`SVD kappa Cancer`),
    ci_lows = c(res$SVD_CI_Normal_Low, res$SVD_CI_Cancer_Low),
    ci_highs = c(res$SVD_CI_Normal_High, res$SVD_CI_Cancer_High),
    group_labels = c("Normal", "Cancer"),
    title = paste("SVD Kappa with 95% CI:", tissue),
    ylab = "SVD κ"
  )
}

# 2. Plot permutation histogram --------------------------------

plot_permutation <- function(tissue, res) {
  make_permutation_plot(
    null_dist = res$perm_dist,
    observed_value = res$`SVD kappa Delta`,
    title = paste("Permutation Test –", tissue),
    xlab = "Δ SVD κ (Cancer - Normal)"
  )
}

# 3. Violin plot for sample-level kappa ------------------------

plot_sample_kappas <- function(tissue, res) {
  values <- c(res$per_sample_n, res$per_sample_c)
  groups <- rep(c("Normal", "Cancer"),
                times = c(length(res$per_sample_n), length(res$per_sample_c)))
  
  make_violin_jitter_plot(
    values = values,
    groups = groups,
    title = paste("Per-Sample SVD Kappa:", tissue),
    ylab = "SVD κ"
  )
}
