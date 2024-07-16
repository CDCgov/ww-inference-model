#' Generates the correlation matrix for a non-spatial process.
#'
#' @param corr_function_params NULL
#'
#' @return Correlation matrix of diagonal ones.
#'
independence_corr_func <- function(
    corr_function_params = list(
      num_sites = NULL
    )) {
  # error managing
  stopifnot(
    "corr.function.params : NOT STRUCTURED CORRECTLY!!!
             List names should be : {num_sites}" =
      names(corr_function_params) == c("num_sites")
  )
  stopifnot(
    "corr.function.params : DOES NOT CONTAIN VALID VALUES!!!
             List should be {int}" =
      is.numeric(corr_function_params$num_sites)
  )

  return(diag(x = 1, nrow = corr_function_params$num_sites))
}
