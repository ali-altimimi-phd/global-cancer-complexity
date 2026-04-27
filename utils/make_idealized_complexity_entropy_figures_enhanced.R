# ------------------------------------------------------------------------------
# File: make_idealized_complexity_entropy_figures_enhanced.R
# Purpose: Generate enhanced idealized 2D schematic plots for latent geometry,
#          complexity, and entropy/chaos regimes.
# Role: Representative figures for write-up; not real data.
# Pipeline: None
# Project: Chaos and Complexity in Cancer
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(MASS)
})

set.seed(123)

# ---- Paths -------------------------------------------------------------------
# This intentionally mirrors the original script's relative project pathing.
# Run from the project root.
OUT_DIR <- file.path("quarto", "resources", "plots", "idealized_2d_enhanced")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Global style -------------------------------------------------------------
COLOR_MAP <- c(
  "Normal" = "#1f77b4",  # blue
  "Tumor"  = "#d62728"   # red
)

AXIS_LIMITS <- c(-4.5, 4.5)
POINT_SIZE  <- 0.5
POINT_ALPHA <- 0.25

# ---- Data builders ------------------------------------------------------------
make_cloud <- function(n = 320, mean = c(0, 0), cov = diag(2), group = "Normal") {
  MASS::mvrnorm(n, mu = mean, Sigma = cov) |>
    as_tibble(.name_repair = "minimal") |>
    setNames(c("z1", "z2")) |>
    mutate(group = group)
}

make_entropy_gain_cloud <- function(n = 320) {
  # Tumor as a noisy, less predictable fan around a transition path.
  t <- runif(n, 0, 1)
  theta <- runif(n, -pi, pi)
  radial_noise <- rnorm(n, mean = 0.25 + 1.35 * t, sd = 0.45 + 0.35 * t)

  tibble(
    z1 = -0.75 + 1.8 * t + radial_noise * cos(theta),
    z2 = -0.45 + 1.4 * t + radial_noise * sin(theta),
    group = "Tumor"
  )
}

make_entropy_loss_cloud <- function(n = 320) {
  # Tumor as constrained/ordered points around a narrow latent trajectory.
  t <- runif(n, 0, 1)
  center_x <- -2.1 + 4.2 * t
  center_y <- -1.6 + 3.2 * t
  tangent <- c(4.2, 3.2) / sqrt(4.2^2 + 3.2^2)
  normal_vec <- c(-tangent[2], tangent[1])
  longitudinal <- rnorm(n, 0, 0.18)
  transverse   <- rnorm(n, 0, 0.16)

  tibble(
    z1 = center_x + longitudinal * tangent[1] + transverse * normal_vec[1],
    z2 = center_y + longitudinal * tangent[2] + transverse * normal_vec[2],
    group = "Tumor"
  )
}

normalize_group_labels <- function(df) {
  df |>
    mutate(group = case_when(
      grepl("Normal", group, ignore.case = TRUE) ~ "Normal",
      TRUE ~ "Tumor"
    ))
}

# ---- Geometry helpers ---------------------------------------------------------
centroid_df <- function(df) {
  df |>
    group_by(group) |>
    summarize(z1 = mean(z1), z2 = mean(z2), .groups = "drop")
}

principal_axes_df <- function(df, scale = 1.7) {
  df |>
    group_by(group) |>
    group_modify(function(.x, .y) {
      cov_mat <- cov(.x[, c("z1", "z2")])
      eig <- eigen(cov_mat)
      vals <- pmax(eig$values, 0)
      vecs <- eig$vectors
      ctr <- colMeans(.x[, c("z1", "z2")])

      map_dfr(seq_along(vals), function(i) {
        v <- vecs[, i]
        len <- scale * sqrt(vals[i])
        tibble(
          axis = paste0("PC", i),
          x = ctr[1] - len * v[1],
          y = ctr[2] - len * v[2],
          xend = ctr[1] + len * v[1],
          yend = ctr[2] + len * v[2]
        )
      })
    }) |>
    ungroup()
}

transition_arrow_df <- function(df) {
  ctr <- centroid_df(df) |>
    pivot_wider(names_from = group, values_from = c(z1, z2))

  tibble(
    x = ctr$z1_Normal,
    y = ctr$z2_Normal,
    xend = ctr$z1_Tumor,
    yend = ctr$z2_Tumor
  )
}

base_plot <- function(df, title, subtitle) {
  df <- normalize_group_labels(df)

  ggplot(df, aes(z1, z2, color = group)) +
    geom_point(alpha = POINT_ALPHA, size = POINT_SIZE) +
    scale_color_manual(values = COLOR_MAP) +
    coord_equal(xlim = AXIS_LIMITS, ylim = AXIS_LIMITS, expand = FALSE) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10.5),
      legend.position = "right",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey88")
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Latent dimension 1",
      y = "Latent dimension 2"
    )
}

add_centroids <- function(p, df) {
  p +
    geom_point(
      data = centroid_df(df),
      aes(z1, z2, fill = group),
      inherit.aes = FALSE,
      shape = 21,
      size = 3.2,
      stroke = 0.7,
      color = "black"
    ) +
    scale_fill_manual(values = COLOR_MAP, guide = "none")
}

add_transition_arrow <- function(p, df) {
  p +
    geom_segment(
      data = transition_arrow_df(df),
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.6,
      alpha = 0.65,
      color = "grey20",
      arrow = arrow(length = unit(0.16, "inches"), type = "closed")
    )
}

save_plot <- function(p, filename, width = 6.2, height = 5.0) {
  ggsave(
    filename = file.path(OUT_DIR, filename),
    plot = p,
    width = width,
    height = height,
    dpi = 300
  )

  invisible(p)
}

# ---- Plot types ---------------------------------------------------------------
plot_latent_geometry <- function(df, title, subtitle, filename, label_text) {
  df <- normalize_group_labels(df)

  p <- base_plot(df, title, subtitle) +
    stat_density_2d(linewidth = 0.35, alpha = 0.30, bins = 4) +
    stat_ellipse(type = "norm", level = 0.82, linewidth = 0.85, alpha = 0.95) +
    annotate("text", x = -4.15, y = 4.10, label = label_text,
             hjust = 0, size = 3.2, color = "grey25")

  p <- add_centroids(p, df)
  p <- add_transition_arrow(p, df)

  save_plot(p, filename)
}

plot_complexity <- function(df, title, subtitle, filename, label_text) {
  df <- normalize_group_labels(df)
  axes <- principal_axes_df(df, scale = 1.55)

  p <- base_plot(df, title, subtitle) +
    stat_ellipse(type = "norm", level = 0.82, linewidth = 0.80, alpha = 0.90) +
    geom_segment(
      data = axes,
      aes(x = x, y = y, xend = xend, yend = yend, color = group),
      inherit.aes = FALSE,
      linewidth = 0.75,
      alpha = 0.80,
      arrow = arrow(length = unit(0.10, "inches"), type = "closed")
    ) +
    annotate("text", x = -4.15, y = 4.10, label = label_text,
             hjust = 0, size = 3.2, color = "grey25")

  p <- add_centroids(p, df)
  p <- add_transition_arrow(p, df)

  save_plot(p, filename)
}

plot_entropy_path <- function(df, title, subtitle, filename, label_text, path_type = c("diffuse", "ordered")) {
  df <- normalize_group_labels(df)
  path_type <- match.arg(path_type)

  tumor_path <- if (path_type == "diffuse") {
    tibble(
      x = c(-1.2, -0.4, 0.35, 1.15),
      y = c(-0.7, -0.15, 0.35, 0.95),
      xend = c(-0.4, 0.35, 1.15, 1.95),
      yend = c(-0.15, 0.35, 0.95, 1.35)
    )
  } else {
    tibble(
      x = c(-2.25, -1.20, -0.10, 1.05),
      y = c(-1.75, -0.95, -0.10, 0.78),
      xend = c(-1.20, -0.10, 1.05, 2.10),
      yend = c(-0.95, -0.10, 0.78, 1.55)
    )
  }

  p <- base_plot(df, title, subtitle) +
    stat_density_2d(linewidth = 0.35, alpha = 0.28, bins = 5) +
    geom_segment(
      data = tumor_path,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.65,
      alpha = 0.75,
      color = COLOR_MAP[["Tumor"]],
      arrow = arrow(length = unit(0.11, "inches"), type = "closed")
    ) +
    annotate("text", x = -4.15, y = 4.10, label = label_text,
             hjust = 0, size = 3.2, color = "grey25")

  p <- add_centroids(p, df)
  p <- add_transition_arrow(p, df)

  save_plot(p, filename)
}

# ---- Schematic datasets -------------------------------------------------------
normal <- make_cloud(
  cov = matrix(c(1.0, 0.0,
                 0.0, 1.0), nrow = 2),
  group = "Normal"
)

complexity_gain <- bind_rows(
  normal,
  make_cloud(
    cov = matrix(c(1.8, 0.15,
                   0.15, 1.55), nrow = 2),
    group = "Tumor"
  )
)

complexity_loss <- bind_rows(
  normal,
  make_cloud(
    cov = matrix(c(2.15, 0.03,
                   0.03, 0.055), nrow = 2),
    group = "Tumor"
  )
)

entropy_gain <- bind_rows(
  normal,
  make_entropy_gain_cloud()
)

entropy_loss <- bind_rows(
  normal,
  make_entropy_loss_cloud()
)

latent_expansion <- bind_rows(
  normal,
  make_cloud(
    cov = matrix(c(2.85, 0.0,
                   0.0, 2.85), nrow = 2),
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

# ---- Generate six separate figures ------------------------------------------
plot_complexity(
  complexity_gain,
  "Idealized complexity gain",
  "Tumor state space recruits additional effective degrees of freedom.",
  "complexity_gain.png",
  "Interpretive cue: effective dimensionality increases"
)

plot_complexity(
  complexity_loss,
  "Idealized complexity loss",
  "Tumor state space collapses toward fewer dominant modes.",
  "complexity_loss.png",
  "Interpretive cue: variance concentrates in fewer axes"
)

plot_entropy_path(
  entropy_gain,
  "Idealized entropy gain / chaotic dispersion",
  "Tumor samples become less predictable around the transition path.",
  "entropy_gain.png",
  "Interpretive cue: diffuse, noisy departure from normal organization",
  path_type = "diffuse"
)

plot_entropy_path(
  entropy_loss,
  "Idealized entropy loss / anti-chaotic constraint",
  "Tumor samples become ordered along a constrained trajectory.",
  "entropy_loss.png",
  "Interpretive cue: constrained, path-like organization",
  path_type = "ordered"
)

plot_latent_geometry(
  latent_expansion,
  "Idealized latent-space expansion",
  "Tumor samples occupy a larger-radius latent region than normal samples.",
  "latent_expansion.png",
  "Interpretive cue: radius/spread increases"
)

plot_latent_geometry(
  latent_contraction,
  "Idealized latent-space contraction",
  "Tumor samples occupy a smaller-radius latent region than normal samples.",
  "latent_contraction.png",
  "Interpretive cue: radius/spread decreases"
)

message("Saved enhanced idealized plots to: ", OUT_DIR)
