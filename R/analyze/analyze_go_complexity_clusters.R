# ------------------------------------------------------------------------------
# File: analyze_go_complexity_clusters.R
# Purpose: Perform semantic similarity clustering of GO terms for a single
#   comparison and summarize complexity-driven cluster structure
# Role: GO complexity clustering utility
# Pipeline: Reporting
# Project: Cancer Complexity Analysis
# Author: Ali M. Al-Timimi
# Created: 2026
# ------------------------------------------------------------------------------

#' Cluster and summarize GO terms for a comparison using complexity metrics
#'
#' Performs semantic similarity clustering of significant GO terms
#' using the Wang method and `GOSemSim`. Summarizes clusters by
#' consensus direction ("gain", "loss", "mixed").
#'
#' @param cmp Character. Comparison code (e.g., `"BR/BRAD"`).
#' @param complexity_df Data frame. Complexity analysis results with GO term data.
#' @param output_dir Character. Directory where the HTML summary table will be written.
#' @param go_mode Character. Optional GO ontology filter: one of `"GO_BP"`, `"GO_MF"`, or `"GO_UNSPECIFIED"`.
#'
#' @return Invisibly returns a data frame summarizing GO clusters.
analyze_go_complexity_clusters <- function(cmp,
                                           complexity_df,
                                           output_dir,
                                           go_mode = NULL) {
  library(GOSemSim)
  library(org.Hs.eg.db)
  
  message("  → Clustering GO complexity terms for: ", cmp,
          if (!is.null(go_mode)) paste0(" (", go_mode, ")") else "")
  
  df <- complexity_df |>
    dplyr::filter(
      comparison == cmp,
      mode == "GO",
      !is.na(gene_set_normalized),
      grepl("^GO:", gene_set_normalized)
    )
  
  if (!is.null(go_mode)) {
    go_mode_short <- stringr::str_replace(go_mode, "GO_", "")
    if (!is.null(go_mode) && go_mode != "GO_UNSPECIFIED") {
      df <- dplyr::filter(df, go_ontology == go_mode_short)
    }
  }
  
  message("  → After ontology filter: ", nrow(df), " GO terms remain.")

  if (nrow(df) < 2) {
    writeLines(
      "<p><em>Not enough GO terms for clustering.</em></p>",
      file.path(output_dir, paste0(sanitize_comparison(cmp), "_go_clusters.html"))
    )
    return(NULL)
  }
  
  go_terms <- unique(df$gene_set_normalized)
  
  ont_short <- dplyr::case_when(
    go_mode == "GO_BP" ~ "BP",
    go_mode == "GO_MF" ~ "MF",
    go_mode == "GO_CC" ~ "CC",
    TRUE ~ "MF"  # default fallback
  )
  
  hsGO <- godata(annoDb = "org.Hs.eg.db", ont = ont_short)
  valid_terms <- intersect(go_terms, names(hsGO@IC))
  
  if (length(valid_terms) < 2)
    return(NULL)
  
  sim_matrix <- mgoSim(
    valid_terms,
    valid_terms,
    semData = hsGO,
    measure = "Wang",
    combine = NULL
  )
  
  sim_matrix[is.na(sim_matrix)] <- 0
  dist_matrix <- as.dist(1 - sim_matrix)
  hc <- hclust(dist_matrix, method = "ward.D2")
  clusters <- cutree(hc, h = 1 - 0.7)
  
  cluster_df <- tibble::tibble(gene_set_normalized = names(clusters),
                               cluster_id = clusters)
  joined <- dplyr::inner_join(df, cluster_df, by = "gene_set_normalized")
  
  cluster_save_dir <- "output/global_cancer/RData/go_clusters"
  fs::dir_create(cluster_save_dir)
  
  save_path <- file.path(cluster_save_dir,
                         paste0(sanitize_comparison(cmp), "_go_clusters.RData"))
  save(joined, file = save_path)
  
  summary_df <- joined |>
    dplyr::group_by(cluster_id) |>
    dplyr::summarise(
      n_terms = dplyr::n(),
      gained = sum(direction == "gained", na.rm = TRUE),
      lost = sum(direction == "lost", na.rm = TRUE),
      direction_consensus = dplyr::case_when(
        gained > lost ~ "gain",
        lost > gained ~ "loss",
        TRUE ~ "mixed"
      ),
      combined_p = if (dplyr::n() == 1) {
        p_perm[1]
      } else {
        valid_p <- p_perm[!is.na(p_perm) & p_perm > 0 & p_perm <= 1]
        if (length(valid_p) == 0) {
          NA_real_
        } else if (length(valid_p) == 1) {
          valid_p[1]
        } else {
          tryCatch({
            metap::sumlog(valid_p)$p
          }, error = function(e) {
            message("sumlog failed: ", e$message)
            NA_real_
          })
        }
      },
      representative_terms = paste(head(gene_set_name, 3), collapse = "; ")
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!(n_terms == 1 & combined_p > 0.05)) |>
    dplyr::filter(!(combined_p > 0.05 & direction_consensus == "mixed")) |>
    dplyr::arrange(combined_p)
  
  html_path <- file.path(output_dir,
                         paste0(sanitize_comparison(cmp), "_go_clusters.html"))
  
  if (nrow(summary_df) == 0) {
    writeLines("<p><em>No significant GO clusters found.</em></p>", html_path)
  } else {
    summary_split <- dplyr::group_split(summary_df, direction_consensus)
    
    html_blocks <- purrr::map_chr(summary_split, function(subgroup) {
      subgroup_title <- stringr::str_to_title(unique(subgroup$direction_consensus))
      
      tbl <- subgroup |>
        dplyr::arrange(combined_p) |>
        dplyr::mutate(`p-value` = sprintf("%.3f", combined_p)) |>
        dplyr::select(
          Cluster = cluster_id,
          `# Terms` = n_terms,
          `p-value`,
          `Representative Terms` = representative_terms
        )
      
      tbl |>
        kableExtra::kable(format = "html", align = "lrrl") |>
        kableExtra::kable_styling(
          bootstrap_options = c("striped", "hover", "condensed"),
          full_width = TRUE
        ) |>
        kableExtra::add_footnote(
          label = "Combined p-values are computed using Fisher’s method (sumlog) across the GO terms in each cluster.",
          notation = "symbol"
        ) |>
        as.character() |>
        (\(html_table) paste0("<h4>", subgroup_title, " Complexity Clusters</h4>\n", html_table))()
    })
    
    writeLines(html_blocks, html_path)
  }
  
  invisible(summary_df)
}
