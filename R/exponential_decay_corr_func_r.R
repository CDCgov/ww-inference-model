#' Generates the correlation matrix for an exponential
#' decay correlation function.
#'
#' @param corr_function_params Structured as follows : dist_matrix, phi, l.
#' Ensure that this list follows this structure. Distance matrix is constructed
#' externally, phi is a numeric value, l is also a numeric value.  Defaults are
#' set to NULL.
#' @export
#'
#' @return Correlation matrix.
#'
exponential_decay_corr_func_r <- function(
    corr_function_params = list(
      dist_matrix = NULL,
      phi = NULL,
      l = NULL
    )) {
  # error managing
  stopifnot(
    "corr_function_params : NOT STRUCTURED CORRECTLY!!!\n
             List names should be : {dist_matrix, phi, l}" =
      names(corr_function_params) == c("dist_matrix", "phi", "l")
  )
  # presets
  dist_matrix <- corr_function_params$dist_matrix
  phi <- corr_function_params$phi
  l <- corr_function_params$l
  stopifnot(
    "corr_function_params : DOES NOT CONTAIN VALID VALUES!!!\n
             List should be {matrix, int, int}" =
      is.matrix(dist_matrix) & is.numeric(phi) & is.numeric(l)
  )
  stopifnot(
    "corr_function_params : DIST.MATRIX NOT SQUARE!!!\n
             Make sure dist_matrix is square" =
      nrow(dist_matrix) == ncol(dist_matrix)
  )


  return(exp(-(dist_matrix / phi)^l))
}
