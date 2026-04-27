library(tidyverse)
library(ggplot2)

simulate_lorenz <- function(
    sigma = 10,
    rho = 28,
    beta = 8/3,
    dt = 0.01,
    n = 6000,
    init = c(x = 1, y = 1, z = 1)
) {
  out <- matrix(NA_real_, nrow = n, ncol = 3)
  colnames(out) <- c("x", "y", "z")
  out[1, ] <- init
  
  for (i in 2:n) {
    x <- out[i - 1, "x"]
    y <- out[i - 1, "y"]
    z <- out[i - 1, "z"]
    
    dx <- sigma * (y - x)
    dy <- x * (rho - z) - y
    dz <- x * y - beta * z
    
    out[i, ] <- c(
      x + dx * dt,
      y + dy * dt,
      z + dz * dt
    )
  }
  
  as_tibble(out) |>
    mutate(step = row_number())
}

plot_attractor_2d <- function(df, title, subtitle = NULL) {
  ggplot(df, aes(x, z)) +
    geom_path(alpha = 0.55, linewidth = 0.25) +
    coord_equal() +
    theme_minimal(base_size = 13) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Idealized state dimension 1",
      y = "Idealized state dimension 2"
    )
}

lorenz <- simulate_lorenz()

chaos_gain <- lorenz |>
  slice(1000:n())

anti_chaos_gain <- lorenz |>
  slice(1000:n()) |>
  mutate(
    x = x * 0.35,
    z = z * 0.35
  )

chaos_loss <- lorenz |>
  slice(1000:n()) |>
  mutate(
    x = x * 0.45,
    z = z * 0.45
  )

anti_chaos_loss <- lorenz |>
  slice(1000:n()) |>
  mutate(
    x = x + rnorm(n(), sd = 6),
    z = z + rnorm(n(), sd = 6)
  )

p_chaos_gain <- plot_attractor_2d(
  chaos_gain,
  "Attractor-inspired schematic: chaos gain",
  "Expanded irregular trajectories suggest increased instability or exploratory state-space behavior."
)

p_anti_chaos_gain <- plot_attractor_2d(
  anti_chaos_gain,
  "Attractor-inspired schematic: anti-chaos gain",
  "Contracted trajectories suggest increased constraint, canalization, or predictability."
)

p_chaos_loss <- plot_attractor_2d(
  chaos_loss,
  "Attractor-inspired schematic: chaos loss",
  "Reduced trajectory extent suggests diminished irregular exploration of state space."
)

p_anti_chaos_loss <- plot_attractor_2d(
  anti_chaos_loss,
  "Attractor-inspired schematic: anti-chaos loss",
  "Noisy dispersion away from the attractor suggests loss of constraint or canalization."
)

dir.create("output/global_cancer/plots/idealized_attractors", recursive = TRUE, showWarnings = FALSE)

ggsave("output/global_cancer/plots/idealized_attractors/chaos_gain_schematic.png", p_chaos_gain, width = 6, height = 5, dpi = 300)
ggsave("output/global_cancer/plots/idealized_attractors/anti_chaos_gain_schematic.png", p_anti_chaos_gain, width = 6, height = 5, dpi = 300)
ggsave("output/global_cancer/plots/idealized_attractors/chaos_loss_schematic.png", p_chaos_loss, width = 6, height = 5, dpi = 300)
ggsave("output/global_cancer/plots/idealized_attractors/anti_chaos_loss_schematic.png", p_anti_chaos_loss, width = 6, height = 5, dpi = 300)