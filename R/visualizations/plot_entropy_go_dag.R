# Status: Legacy / exploratory visualization
# Note:
# This script predates the unified GO semantic summarization framework and
# operates on deprecated entropy cluster outputs.

#' Build GO DAG for clustered entropy terms and visualize it
#'
#' @param cmp Character. Comparison code (e.g., "PB/T-ALL") or category name if `category_mode = TRUE`.
#' @param ontology Character. GO ontology ("MF", "BP", "CC"). Default is "MF".
#' @param cluster_dir Character. Path to .RData files containing clustered entropy results.
#' @param palette Named vector mapping spectral directions to colors.
#' @param category_mode Logical. If TRUE, load and merge all comparisons for a tumor category.
#'
#' @return A `ggraph` plot object or NULL if DAG cannot be built.

plot_entropy_go_dag <- function(
    cmp,
    ontology = "BP",
    cluster_dir = "output/global_cancer/RData/go_entropy_clusters",
    palette = c(
      "Strongly Chaotic" = "#d73027",
      "Moderately Chaotic" = "#fc8d59",
      "Neutral" = "#999999",
      "Moderately Ordered" = "#91bfdb",
      "Strongly Ordered" = "#4575b4"
    ),
    category_mode = FALSE
) {
  library(GO.db)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  
  # Helper to build DAG from GO.db
  build_go_dag <- function(go_ids, ontology = "MF") {
    go_anno <- AnnotationDbi::select(
      GO.db,
      keys = go_ids,
      columns = c("GOID", "ONTOLOGY", "TERM"),
      keytype = "GOID"
    )
    
    filtered_ids <- go_anno |>
      dplyr::filter(ONTOLOGY == ontology) |>
      dplyr::pull(GOID) |>
      unique()
    
    if (length(filtered_ids) < 2) return(NULL)
    
    # Collect parent-child edges
    all_edges <- purrr::map_dfr(filtered_ids, function(goid) {
      parents <- tryCatch(AnnotationDbi::get(goid, GO.db::GOPARENTS), error = function(e) NA)
      if (is.na(parents)) return(NULL)
      tibble(from = parents, to = goid)
    })
    
    if (nrow(all_edges) == 0) return(NULL)
    
    all_nodes <- unique(c(all_edges$from, all_edges$to))
    g <- igraph::graph_from_data_frame(all_edges, directed = TRUE, vertices = all_nodes)
    return(g)
  }
  
  # Load data: one or more .RData files depending on mode
  load_joined <- function(file) {
    load(file)
    if (!exists("joined")) stop("Missing 'joined' object in file: ", file)
    joined
  }
  
  if (category_mode) {
    # Load all comparisons for category
    files <- list.files(cluster_dir, pattern = "\\.RData$", full.names = TRUE)
    files <- files[grepl(sanitize_comparison(cmp), files)]
    if (length(files) == 0) stop("No matching files for category: ", cmp)
    joined_all <- purrr::map_dfr(files, load_joined)
  } else {
    fname <- file.path(cluster_dir, paste0(sanitize_comparison(cmp), "_go_entropy_clusters.RData"))
    if (!file.exists(fname)) stop("File not found: ", fname)
    joined_all <- load_joined(fname)
  }
  
  # Extract valid GO terms
  go_ids <- unique(joined_all$gene_set_normalized)
  ontology_info <- AnnotationDbi::select(
    GO.db,
    keys = go_ids,
    columns = c("ONTOLOGY", "TERM"),
    keytype = "GOID"
  )
  
  filtered_ids <- ontology_info |>
    dplyr::filter(ONTOLOGY == ontology, !is.na(TERM)) |>
    dplyr::pull(GOID) |>
    unique()
  
  if (length(filtered_ids) < 2) {
    message("Not enough valid GO terms in ontology '", ontology, "'")
    return(invisible(NULL))
  }
  
  # Build graph
  gograph <- build_go_dag(filtered_ids, ontology)
  if (is.null(gograph)) {
    message("No parent-child edges available among GO terms.")
    return(invisible(NULL))
  }
  
  # Annotate direction for each node
  direction_df <- joined_all |>
    dplyr::filter(gene_set_normalized %in% filtered_ids) |>
    dplyr::distinct(gene_set_normalized, spectral_direction)
  
  V(gograph)$direction <- direction_df$spectral_direction[match(names(V(gograph)), direction_df$gene_set_normalized)]
  V(gograph)$direction[is.na(V(gograph)$direction)] <- "Neutral"
  
  # Plot
  title_str <- if (category_mode) paste("GO DAG for category:", cmp) else paste("GO DAG for", cmp)
  
  ggraph(gograph, layout = "dendrogram", circular = FALSE) +
    geom_edge_elbow(color = "grey80") +
    geom_node_point(aes(color = direction), size = 3) +
    geom_node_text(aes(label = name), hjust = 1.1, size = 3, vjust = 0.3, check_overlap = TRUE) +
    scale_color_manual(values = palette, name = "Entropy Direction") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(title = title_str, x = NULL, y = NULL)
}
