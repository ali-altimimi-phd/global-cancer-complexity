run_latent_notebook <- function(notebook_name, logger) {
  notebook_path <- file.path(latent_notebook_dir, notebook_name)
  
  if (!file.exists(notebook_path)) {
    logger$log(sprintf("Notebook not found: %s", notebook_path), section = "ERROR")
    stop(sprintf("Notebook not found: %s", notebook_path), call. = FALSE)
  }
  
  if (isTRUE(execute_notebooks_inplace)) {
    output_args <- c("--inplace")
    output_label <- notebook_path
  } else {
    output_path <- file.path(notebook_executed_dir, notebook_name)
    output_args <- c("--output", shQuote(output_path))
    output_label <- output_path
  }
  
  args <- c(
    jupyter_prefix_args,
    "nbconvert",
    "--to", "notebook",
    "--execute",
    sprintf("--ExecutePreprocessor.timeout=%s", notebook_timeout_seconds),
    output_args,
    shQuote(notebook_path)
  )
  
  logger$timed(sprintf("Notebook %s", notebook_name), {
    logger$log(sprintf("Launching notebook: %s", notebook_path), section = "NOTEBOOKS")
    
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(latent_notebook_dir)
    
    output <- tryCatch(
      system2(command = jupyter_command, args = args, stdout = TRUE, stderr = TRUE),
      warning = function(w) structure(conditionMessage(w), status = 1L),
      error = function(e) structure(conditionMessage(e), status = 1L)
    )
    
    if (length(output) > 0) {
      for (line in output) logger$log(as.character(line), section = "JUPYTER")
    }
    
    exit_status <- attr(output, "status")
    if (is.null(exit_status)) exit_status <- 0L
    
    if (!identical(as.integer(exit_status), 0L)) {
      logger$log(sprintf("Notebook failed with exit status %s: %s", exit_status, notebook_name), section = "ERROR")
      stop(sprintf("Notebook failed: %s", notebook_name), call. = FALSE)
    }
    
    logger$log(sprintf("Completed notebook successfully: %s", output_label), section = "NOTEBOOKS")
  })
  
  invisible(TRUE)
}


run_latent_notebooks_stage <- function(logger) {
  logger$log("📓 Executing latent-space analysis notebooks...", section = "NOTEBOOKS")
  
  for (notebook_name in latent_notebooks) {
    run_latent_notebook(notebook_name, logger = logger)
  }
  
  logger$log("✅ Latent-space notebooks completed successfully.", section = "NOTEBOOKS")
}