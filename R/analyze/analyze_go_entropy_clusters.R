# ------------------------------------------------------------------------------
# File: analyze_go_entropy_clusters.R
# Purpose: Perform semantic similarity clustering of GO terms for a single
#   comparison and summarize entropy-driven cluster structure
# Role: GO entropy clustering utility
# Pipeline: Reporting
# Project: Cancer Complexity Analysis
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------

# ---- STATUS ------------------------------------------------------------------
# Status: Inactive; pending validation

#' Cluster and summarize GO terms for a comparison using entropy metrics
#'
#' This function is intended to parallel `analyze_go_complexity_clusters()` but
#' is not currently invoked by the active reporting workflow.
#'
#' Retained for development and validation.
analyze_go_entropy_clusters <- function(cmp, entropy_df, output_dir = "quarto/resources/tables/go_entropy") {
  library(GOSemSim)
  library(org.Hs.eg.db)
  
  fs::dir_create(output_dir)
  
  df <- entropy_df |>
    dplyr::filter(
      comparison == cmp,
      mode == "GO",
      !is.na(gene_set_normalized),
      grepl("^GO:", gene_set_normalized),
      p_perm <= 0.05
    )
  
  if (nrow(df) < 2) {
    writeLines(
      "<p><em>Not enough significant GO entropy terms for clustering.</em></p>",
      file.path(output_dir, paste0(sanitize_comparison(cmp), "_go_entropy_clusters.html"))
    )
    return(NULL)
  }
  
  go_terms <- unique(df$gene_set_normalized)
  hsGO <- godata(annoDb = "org.Hs.eg.db", ont = "MF")
  valid_terms <- intersect(go_terms, names(hsGO@IC))
  if (length(valid_terms) < 2) return(NULL)
  
  sim_matrix <- mgoSim(valid_terms, valid_terms, semData = hsGO, measure = "Wang", combine = NULL)
  sim_matrix[is.na(sim_matrix)] <- 0
  dist_matrix <- as.dist(1 - sim_matrix)
  hc <- hclust(dist_matrix, method = "ward.D2")
  clusters <- cutree(hc, h = 1 - 0.7)
  
  cluster_df <- tibble::tibble(gene_set_normalized = names(clusters), cluster_id = clusters)
  joined <- dplyr::inner_join(df, cluster_df, by = "gene_set_normalized")
  
  # ---- Save full joined cluster data to RData file ----
  cluster_save_dir <- "output/global_cancer/RData/go_entropy_clusters"
  fs::dir_create(cluster_save_dir)
  
  save_path <- file.path(cluster_save_dir, paste0(sanitize_comparison(cmp), "_go_entropy_clusters.RData"))
  save(joined, file = save_path)
  # -----------------------------------------------------
  
  summary_df <- joined |>
    dplyr::group_by(cluster_id) |>
    dplyr::summarise(
      n_terms = dplyr::n(),
      combined_p = {
        valid_p <- p_perm[!is.na(p_perm) & p_perm > 0 & p_perm <= 1]
        if (length(valid_p) == 0) {
          NA_real_
        } else if (length(valid_p) == 1) {
          valid_p[1]
        } else {
          tryCatch({
            res <- metap::sumlog(valid_p)$p
            message("Entropy cluster ", dplyr::cur_group_id(), ": combined p = ", round(res,3))
            res
          }, error = function(e) {
            message("sumlog failed: ", e$message)
            NA_real_
          })
        }
      },
      direction_consensus = {
        cats <- tolower(spectral_direction)
        top <- names(which.max(table(cats)))
        stringr::str_to_title(top)
      },
      representative_terms = paste(head(gene_set_name, 3), collapse = "; ")
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!(n_terms == 1 & combined_p > 0.05)) |>
    dplyr::filter(!(combined_p > 0.05 & direction_consensus == "Neutral")) |>
    dplyr::arrange(combined_p)
  
  html_path <- file.path(output_dir, paste0(sanitize_comparison(cmp), "_go_entropy_clusters.html"))
  
  if (nrow(summary_df) == 0) {
    writeLines("<p><em>No significant GO entropy clusters found.</em></p>", html_path)
  } else {
    summary_split <- dplyr::group_split(summary_df, direction_consensus)
    
    html_blocks <- purrr::map_chr(summary_split, function(subgroup) {
      subgroup_title <- unique(subgroup$direction_consensus)
      
      tbl <- subgroup |>
        dplyr::arrange(combined_p) |>
        dplyr::mutate(`p-value` = sprintf("%.3f", combined_p)) |>
        dplyr::select(
          Cluster = cluster_id,
          `# Terms` = n_terms,
          `p-value`,
          `Representative Terms` = representative_terms
        )
      
      html_table <- tbl |>
        kableExtra::kable("html", align = "lrrl") |>
        kableExtra::kable_styling(
          bootstrap_options = c("striped", "hover", "condensed"),
          full_width = TRUE
        ) |>
        kableExtra::add_footnote(
          label = "Combined p-values are computed using Fisher’s method (sumlog) across the GO terms in each entropy cluster.",
          notation = "symbol"
        ) |>
        as.character()
      
      paste0("<h4>", subgroup_title, " Entropy Clusters</h4>\n", html_table)
    })
    
    writeLines(html_blocks, html_path)
  }
  
  invisible(summary_df)
}
