

save_all_complexity_plots <- function(df, folder) {
  for (row in seq_len(nrow(df))) {
    res <- df$tidy_results[[row]]
    tissue <- df$Tissue[[row]]
    
    plots <- preview_complexity_plots(res, tissue)
    
    ggsave(here::here("output", "plots", folder, glue::glue("kappa_ci_{tissue}.png")),
           plot = plots$kappa_ci, width = 6, height = 4, bg = "white")
    
    ggsave(here::here("output", "plots", folder, glue::glue("perm_test_{tissue}.png")),
           plot = plots$permutation, width = 6, height = 4, bg = "white")
    
    ggsave(here::here("output", "plots", folder, glue::glue("sample_kappas_{tissue}.png")),
           plot = plots$sample_kappas, width = 6, height = 4, bg = "white")
  }
}

# Call function for carcinomas, etc.
# save_all_complexity_plots(df_complexity_carcinomas, "carcinomas")

