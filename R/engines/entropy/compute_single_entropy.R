# -------------------------------------------------------------
# Compute entropy metrics for individual matrices
# -------------------------------------------------------------

compute_matrix_entropy <- function(label, chip, mat) {
  stopifnot(is.matrix(mat))
  
  tibble::tibble(
    Label    = label,
    Chip     = chip,
    Samples  = ncol(mat),
    Shannon  = compute_shannon_entropy(mat),
    Spectral = compute_spectral_entropy(mat)
  )
}
