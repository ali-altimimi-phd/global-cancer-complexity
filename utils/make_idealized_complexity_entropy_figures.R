# ------------------------------------------------------------------------------
# File: make_idealized_complexity_entropy_plots.R
# Purpose: Generate idealized 2D schematic plots for complexity/entropy regimes
# Role: Ued in write-up
# Pipeline: None
# Project: Chaos and Complexity in Cancer
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------

library(tidyverse)
library(MASS)

set.seed(123)

OUT_DIR <- file.path("quarto", "resources", "plots", "idealized_2d")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

COLOR_MAP <- c(
  "Normal" = "#1f77b4",  # blue
  "Tumor"  = "#d62728"   # red
)

make_cloud <- function(n = 300, mean = c(0, 0), cov = diag(2), group = "Normal") {
  MASS::mvrnorm(n, mu = mean, Sigma = cov) |>
    as_tibble(.name_repair = "minimal") |>
    setNames(c("z1", "z2")) |>
    mutate(group = group)
}

normalize_group_labels <- function(df) {
  df |>
    mutate(group = case_when(
      grepl("Normal", group, ignore.case = TRUE) ~ "Normal",
      TRUE ~ "Tumor"
    ))
}

plot_and_save <- function(df, title, subtitle, filename) {
  df <- normalize_group_labels(df)
  
  p <- ggplot(df, aes(z1, z2, color = group)) +
    geom_point(alpha = 0.60, size = 2) +
    scale_color_manual(values = COLOR_MAP) +
    coord_equal() +
    theme_minimal(base_size = 14) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Latent dimension 1",
      y = "Latent dimension 2",
      color = NULL
    )
  
  ggsave(
    filename = file.path(OUT_DIR, filename),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
  
  return(p)
}

normal <- make_cloud(
  cov = matrix(c(1.0, 0.0,
                 0.0, 1.0), nrow = 2),
  group = "Normal"
)

complexity_gain <- bind_rows(
  normal,
  make_cloud(
    cov = matrix(c(2.2, 0.6,
                   0.6, 1.8), nrow = 2),
    group = "Tumor"
  )
)

complexity_loss <- bind_rows(
  normal,
  make_cloud(
    cov = matrix(c(2.2, 0.05,
                   0.05, 0.05), nrow = 2),
    group = "Tumor"
  )
)

entropy_gain <- bind_rows(
  normal,
  make_cloud(
    cov = matrix(c(2.5, 0.0,
                   0.0, 2.5), nrow = 2),
    group = "Tumor"
  )
)

entropy_loss <- bind_rows(
  normal,
  make_cloud(
    cov = matrix(c(0.25, 0.0,
                   0.0, 0.25), nrow = 2),
    group = "Tumor"
  )
)

latent_expansion <- bind_rows(
  normal,
  make_cloud(
    cov = matrix(c(3.0, 0.0,
                   0.0, 3.0), nrow = 2),
    group = "Tumor"
  )
)

latent_contraction <- bind_rows(
  normal,
  make_cloud(
    cov = matrix(c(0.20, 0.0,
                   0.0, 0.20), nrow = 2),
    group = "Tumor"
  )
)

p_complexity_gain <- plot_and_save(
  complexity_gain,
  "Idealized complexity gain",
  "Tumor samples occupy a broader, higher-dimensional state space.",
  "complexity_gain.png"
)

p_complexity_loss <- plot_and_save(
  complexity_loss,
  "Idealized complexity loss",
  "Tumor samples collapse toward a lower-dimensional structure.",
  "complexity_loss.png"
)

p_entropy_gain <- plot_and_save(
  entropy_gain,
  "Idealized entropy gain",
  "Tumor samples become more diffuse and less predictable.",
  "entropy_gain.png"
)

p_entropy_loss <- plot_and_save(
  entropy_loss,
  "Idealized entropy loss",
  "Tumor samples become more concentrated and ordered.",
  "entropy_loss.png"
)

p_latent_expansion <- plot_and_save(
  latent_expansion,
  "Idealized latent-space expansion",
  "Tumor samples occupy a larger-radius latent region than normal samples.",
  "latent_expansion.png"
)

p_latent_contraction <- plot_and_save(
  latent_contraction,
  "Idealized latent-space contraction",
  "Tumor samples occupy a smaller-radius latent region than normal samples.",
  "latent_contraction.png"
)

message("Saved idealized plots to: ", OUT_DIR)