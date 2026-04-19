library(ggplot2)
library(dplyr)
library(tidyr)
source(here::here("R/helpers/category_utils.R"))
source(here::here("R/helpers/plot_utils.R"))  # for save_png_light()

#' Generate and save heatmaps of directionality and confidence for each metric
#'
#' @param df Summary data frame with weighted scores and confidence labels
#' @param categories Vector of cancer categories to include (e.g., c("carcinomas", "blastomas"))
#' @param output_dir Directory to save the PNGs
#' @export
plot_category_metric_heatmaps <- function(df, categories, output_dir = plots_dir) {
  df <- add_cancer_category(df)
  
  df_long <- df %>%
    filter(cancer_category %in% categories) %>%
    pivot_longer(
      cols = c(complexity_weighted_score, entropy_weighted_score_shannon, entropy_weighted_score_spectral),
      names_to = "metric",
      values_to = "score"
    ) %>%
    mutate(
      confidence = case_when(
        metric == "complexity_weighted_score" ~ complexity_confidence_label,
        metric == "entropy_weighted_score_shannon" ~ confidence_label_shannon,
        metric == "entropy_weighted_score_spectral" ~ confidence_label_spectral
      ),
      metric = recode(metric,
                      complexity_weighted_score = "Complexity",
                      entropy_weighted_score_shannon = "Shannon Entropy",
                      entropy_weighted_score_spectral = "Spectral Entropy")
    )
  
  for (category in categories) {
    df_cat <- df_long %>% filter(cancer_category == category)
    
    p <- ggplot(df_cat, aes(x = comparison, y = mode_label, fill = score)) +
      geom_tile(color = "white") +
      facet_wrap(~ metric, ncol = 1) +
      geom_text(aes(label = confidence), size = 2.5) +
      scale_fill_gradient2(
        low = "blue", high = "red", mid = "white", midpoint = 0,
        name = "Weighted Score"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold")
      ) +
      labs(
        title = paste("Metric Direction & Confidence –", category),
        x = "Comparison", y = "Data Mode"
      )
    
    out_file <- file.path(output_dir, paste0("heatmap_", category, ".png"))
    save_png_light(p, out_file, width = 8, height = 10)
  }
}
