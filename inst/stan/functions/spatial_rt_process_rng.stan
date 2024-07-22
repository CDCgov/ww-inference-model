/**
  * Generates the subpopulation Rt trajectories for a spatially correlated model
  * @param log_state_rt "State" level unadjusted Rt, on log scale.
  * @param phi_rt Coefficient for AR(1) temporal correlation on
  * subpopulation deviations.
  * @param sigma_eps Parameter for construction of covariance matrix
  * of spatial epsilon.
  * @param n_sites Number of sites/subpopulations
  * @param delta_0 Initial value of delta
  * @return A matrix for unadjusted Rt where rows are subpopulations
  * and columns are time points.
  */
matrix spatial_rt_process_rng(vector log_state_rt, real phi_rt, real sigma_eps, int n_sites, vector delta_0) {


  // correlation matrix constr.
  matrix[n_sites,n_sites] omega_matrix_eps = independence_corr_func(n_sites);
  matrix[n_sites,n_sites] sigma_matrix_eps = sigma_eps^2 * omega_matrix_eps;

  // presets
  int n_time = dims(log_state_rt)[1];
  matrix[n_sites,n_time] log_site_rt;
  matrix[n_sites,n_time] delta;
  vector[n_sites] eps_vec;

  delta[,1] = delta_0;
    for (t_i in 2:n_time) {
      eps_vec = multi_normal_rng(rep_vector(0,n_sites), sigma_matrix_eps);
      delta[,t_i] = phi_rt * delta[,t_i - 1] + eps_vec;
    }

  for (t_i in 1:n_time) {
    log_site_rt[,t_i] = (log_state_rt[t_i] * rep_vector(1,n_sites)) + delta[,t_i];
  }

  return log_site_rt;
}
