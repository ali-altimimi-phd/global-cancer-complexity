library(dplyr)
library(ggplot2)
library(tidyr)
library(glue)
library(purrr)
library(knitr)
library(kableExtra)

# Step 1: Extract filtered probe counts
extract_probe_counts_df <- function(chips = c("hu35ksuba", "hu6800")) {
  purrr::map_dfr(chips, function(chip) {
    res_name <- glue("res_{chip}")
    if (!exists(res_name, envir = .GlobalEnv)) return(NULL)
    
    res <- get(res_name, envir = .GlobalEnv)
    
    purrr::imap_dfr(res, function(group_data, group) {
      if (group == "__summary__") return(NULL)
      
      purrr::imap_dfr(group_data, function(cmp, cmp_name) {
        tibble(
          chip = chip,
          group = group,
          comparison = cmp_name,
          n_filtered = length(cmp$filtered_probes)
        )
      })
    })
  })
}

# Step 2: Generate formatted kable table for one group
make_probe_table <- function(df, group_name) {
  df %>%
    filter(group == group_name) %>%
    arrange(comparison, chip) %>%
    select(Chip = chip, Comparison = comparison, `Filtered Probes` = n_filtered) %>%
    kable("html", caption = glue("Filtered Probe Counts for {group_name} Group")) %>%
    kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
}

# Step 3: Generate ggplot bar chart for one group
make_probe_plot <- function(df, group_name) {
  df_group <- df %>% filter(group == group_name)
  
  ggplot(df_group, aes(x = comparison, y = n_filtered, fill = chip)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    labs(
      title = glue("Filtered Probes per Comparison: {group_name}"),
      x = "Comparison",
      y = "Filtered Probes",
      fill = "Chip"
    ) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Optional helper to list available groups
get_probe_groups <- function(chips = c("hu35ksuba", "hu6800")) {
  df <- extract_probe_counts_df(chips)
  unique(df$group)
}

save_probe_plot <- function(df, group_name, out_dir = "output/global_cancer/plots/probes") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  p <- make_probe_plot(df, group_name)
  out_path <- file.path(out_dir, glue::glue("{group_name}_barplot.png"))
  
  ggsave(out_path, plot = p, width = 7, height = 4.5, dpi = 300)
  message(glue("✅ Saved: {out_path}"))
}
