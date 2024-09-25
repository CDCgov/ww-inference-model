#' Passes user specified correlation matrix, used for `generate_simulated_data`.
#'
#' @param corr_function_params Structured as follows : {corr_matrix}.
#' Ensure that this list follows this structure. Correlation matrix with proprer
#' dimensions.  Defaults are set to NULL.
#' @export
#'
#' @return The same correlation matrix that was passed in.
#'
custom_corr_matrix_func <- function(
    corr_function_params = list(
      corr_matrix = NULL
    )) {
  # error managing
  stopifnot(
    "corr_function_params : NOT STRUCTURED CORRECTLY!!!\n
             List names should be : {corr_matrix}" =
      names(corr_function_params) == c("corr_matrix")
  )

  return(corr_function_params$corr_matrix)
}
