#' Generates a random correlation matrix.
#'
#' @param corr_function_params Structured as follows : {n_sites, seed}.
#' Ensure that this list follows this structure. Number of sites n_sites, and
#' seed desired for rng.  Defaults are set to NULL.
#' @export
#'
#' @return Correlation matrix.  A symmetric matrix with off diagonal values
#' sampled randomly from -1 to 1, using a uniform distribution.
#'
rand_corr_matrix_func <- function(
    corr_function_params = list(
      n_sites = NULL,
      seed = NULL
    )) {
  # error managing
  stopifnot(
    "corr_function_params : NOT STRUCTURED CORRECTLY!!!\n
             List names should be : {n_sites, seed}" =
      names(corr_function_params) == c("n_sites", "seed")
  )
  # presets
  n_sites <- corr_function_params$n_sites
  seed <- corr_function_params$seed
  set.seed(seed)
  random_matrix <- matrix(
    data = runif(
      n_sites * n_sites,
      min = -1, max = 1
    ),
    nrow = n_sites
  )
  sym_matrix <- (random_matrix + t(random_matrix)) / 2
  eigen_decomp <- eigen(
    sym_matrix
  )
  positive_semi_definite_matrix <- eigen_decomp$vectors %*%
    diag(
      pmax(
        eigen_decomp$values,
        0
      )
    ) %*%
    t(eigen_decomp$vectors)
  diag(positive_semi_definite_matrix) <- 1


  return(round(positive_semi_definite_matrix, 4))
}
