#' Load Matrix and Comparison Maps into Global Environment
#'
#' Loads objects like `matrices_hu35ksuba`, `matrices_hu6800`,
#' `comparison_map_hu35ksuba`, and `comparison_map_hu6800` from
#' `global_cancer_matrix_maps.RData` into the global environment.
#'
#' @param matrices_path Path to the `.RData` file storing matrix/comparison maps.
#' @param overwrite Logical. If `TRUE`, overwrite existing objects.
#' @export
load_matrix_maps <- function(matrices_path, overwrite = FALSE) {
  stopifnot(file.exists(matrices_path))
  
  loaded_objs <- load(matrices_path)
  for (obj in loaded_objs) {
    if (exists(obj, envir = .GlobalEnv)) {
      if (overwrite) {
        assign(obj, get(obj), envir = .GlobalEnv)
        message(glue::glue("🔁 Overwritten: {obj}"))
      } else {
        message(glue::glue("⏭️ Skipped (already exists): {obj}"))
      }
    } else {
      assign(obj, get(obj), envir = .GlobalEnv)
      message(glue::glue("✅ Loaded: {obj}"))
    }
  }
}

#' Load Filtered Probe Results into Global Environment
#'
#' Loads probe filtering `.rds` results (Stage 3 of analysis pipeline) 
#' for both chips into the global environment as `res_hu35ksuba` and `res_hu6800`.
#' Skips loading if those variables already exist unless `overwrite = TRUE`.
#'
#' @param filtered_dir Path to the directory containing filtered `.rds` files.
#' @param overwrite Logical. If `TRUE`, overwrite existing variables.
load_filtered_results <- function(filtered_dir, overwrite = FALSE) {
  stopifnot(dir.exists(filtered_dir))
  
  files_to_load <- list(
    res_hu35ksuba = file.path(filtered_dir, "filtered_probes_hu35ksuba.rds"),
    res_hu6800 = file.path(filtered_dir, "filtered_probes_hu6800.rds")
  )
  
  for (var_name in names(files_to_load)) {
    file_path <- files_to_load[[var_name]]
    
    if (!file.exists(file_path)) {
      warning(glue::glue("⚠️ File not found: {basename(file_path)}"))
      next
    }
    
    if (exists(var_name, envir = .GlobalEnv)) {
      if (overwrite) {
        assign(var_name, readRDS(file_path), envir = .GlobalEnv)
        message(glue::glue("🔁 Overwritten: {var_name}"))
      } else {
        message(glue::glue("⏭️ Skipped (already exists): {var_name}"))
      }
    } else {
      assign(var_name, readRDS(file_path), envir = .GlobalEnv)
      message(glue::glue("✅ Loaded: {var_name}"))
    }
  }
}

#' Load Comparison Results into Global Environment
#'
#' Loads all `.rds` files from `data_dir` that match the pattern
#' `complexity_results_{something}.rds` or `entropy_results_{something}.rds`.
#' Skips variables that already exist unless `overwrite = TRUE`.
#'
#' @param data_dir Path to the directory containing result `.rds` files.
#' @param overwrite Logical. If `TRUE`, overwrite existing variables in the global environment.
#' @export
load_comparison_results <- function(data_dir, overwrite = FALSE) {
  stopifnot(dir.exists(data_dir))
  
  rds_files <- list.files(
    path = data_dir,
    pattern = "^(complexity|entropy)_results_.*\\.rds$",
    full.names = TRUE
  )
  
  for (file in rds_files) {
    var_name <- tools::file_path_sans_ext(basename(file))
    
    if (exists(var_name, envir = .GlobalEnv)) {
      if (overwrite) {
        assign(var_name, readRDS(file), envir = .GlobalEnv)
        message(glue::glue("🔁 Overwritten: {var_name}"))
      } else {
        message(glue::glue("⏭️ Skipped (already exists): {var_name}"))
      }
    } else {
      assign(var_name, readRDS(file), envir = .GlobalEnv)
      message(glue::glue("✅ Loaded: {var_name}"))
    }
  }
}

#' Load Cleaned Summary Results into Global Environment
#'
#' Loads `complexity_cleaned.rds` and `entropy_cleaned.rds` from `cleaned_dir`
#' into the global environment as `complexity_df` and `entropy_df`, respectively.
#' Skips loading if those variables already exist unless `overwrite = TRUE`.
#'
#' @param cleaned_dir Path to the directory containing cleaned `.rds` files.
#' @param overwrite Logical. If `TRUE`, overwrite existing `complexity_df` and `entropy_df`.
#' @export
load_cleaned_results <- function(cleaned_dir, overwrite = FALSE) {
  stopifnot(dir.exists(cleaned_dir))
  
  files_to_load <- list(
    complexity_df = file.path(aggregate_dir, "complexity_aggregated_results.rds"),
    entropy_df = file.path(aggregate_dir, "entropy_aggregated_results.rds")
  )
  
  for (var_name in names(files_to_load)) {
    file_path <- files_to_load[[var_name]]
    
    if (!file.exists(file_path)) {
      warning(glue::glue("⚠️ File not found: {basename(file_path)}"))
      next
    }
    
    if (exists(var_name, envir = .GlobalEnv)) {
      if (overwrite) {
        assign(var_name, readRDS(file_path), envir = .GlobalEnv)
        message(glue::glue("🔁 Overwritten: {var_name}"))
      } else {
        message(glue::glue("⏭️ Skipped (already exists): {var_name}"))
      }
    } else {
      assign(var_name, readRDS(file_path), envir = .GlobalEnv)
      message(glue::glue("✅ Loaded: {var_name}"))
    }
  }
}

#' Load Combined Summary Results into Global Environment
#'
#' Loads the summary data frame (typically named `summaries_combined_df`)
#' from the configured `.rds` path into the global environment.
#'
#' @param summary_path Path to the `.rds` file storing the combined summary.
#' @param overwrite Logical. If `TRUE`, overwrite `summaries_combined_df` if it exists.
#' @export
load_summaries_df <- function(summaries_dir, overwrite = FALSE) {
  stopifnot(dir.exists(summaries_dir))
  
  var_name <- "summaries_combined_df"
  
  if (exists(var_name, envir = .GlobalEnv)) {
    if (overwrite) {
      assign(var_name, readRDS(summaries_df_path), envir = .GlobalEnv)
      message(glue::glue("🔁 Overwritten: {var_name}"))
    } else {
      message(glue::glue("⏭️ Skipped (already exists): {var_name}"))
    }
  } else {
    assign(var_name, readRDS(summaries_df_path), envir = .GlobalEnv)
    message(glue::glue("✅ Loaded: {var_name}"))
  }
}
