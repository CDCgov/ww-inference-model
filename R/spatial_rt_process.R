#' Generates the subpopulation Rt trajectories for a spatially correlated model
#'
#' @param log_state_rt "State" level unadjusted Rt, on log scale.
#' @param corr_function Function for the correlation structure desired.
#' @param corr_function_params List of parameters of desired correlation
#' function, should contain dist matrix.
#' @param phi_rt Coefficient for AR(1) temporal correlation on
#' subpopulation deviations.
#' @param sigma_eps Parameter for construction of covariance matrix
#' of spatial epsilon.
#' @export
#'
#' @return A matrix for unadjusted Rt where rows are subpopulations
#' and columns are time points.
#'
spatial_rt_process <- function(log_state_rt,
                               corr_function,
                               corr_function_params,
                               phi_rt,
                               sigma_eps) {
  # correlation matrix constr.
  omega_matrix_eps <- corr_function(corr_function_params)
  sigma_matrix_eps <- sigma_eps^2 * omega_matrix_eps

  # presets
  n_subpopulations <- nrow(sigma_matrix_eps)
  n_time <- length(log_state_rt)
  log_site_rt <- matrix(data = 0, nrow = n_subpopulations, ncol = n_time)
  delta <- matrix(data = 0, nrow = n_subpopulations, ncol = n_time)

  # delta constr.
  delta[, 1] <- mvrnorm(
    n = 1,
    mu = matrix(data = 0, nrow = 1, ncol = n_subpopulations),
    Sigma = sigma_matrix_eps
  )
  for (t_i in 2:n_time) {
    eps_vec <- mvrnorm(
      n = 1,
      mu = matrix(data = 0, nrow = 1, ncol = n_subpopulations),
      Sigma = sigma_matrix_eps
    )
    delta[, t_i] <- phi_rt * delta[, t_i - 1] + eps_vec
  }

  # Subpopulation unadjusted Rt constr.
  for (t_i in 1:n_time) {
    log_site_rt[, t_i] <- (log_state_rt[t_i] * matrix(
      data = 1,
      nrow = n_subpopulations,
      ncol = 1
    )) +
      delta[, t_i]
  }

  return(log_site_rt)
}
