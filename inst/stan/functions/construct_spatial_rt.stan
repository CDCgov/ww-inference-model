/**
  * Generates the subpopulation Rt trajectories for a spatially correlated model
  * given all the random components.
  * @param log_state_rt "State" level unadjusted Rt, on log scale.
  * @param spatial_deviation_ar_coeff Coefficient for AR(1) temporal correlation on
  * subpopulation deviations.
  * @param spatial_deviation_noise_matrix matrix with random vectors of noise
  * for spatial devaiation process.
  * @return A matrix for unadjusted Rt where rows are subpopulations
  * and columns are time points.
  */
matrix construct_spatial_rt(vector log_state_rt,
                            real spatial_deviation_ar_coeff,
                            matrix spatial_deviation_noise_matrix) {


  // presets
  int n_time = dims(log_state_rt)[1];
  int n_sites = rows(spatial_deviation_noise_matrix);
  matrix[n_sites,n_time] log_site_rt;
  vector[n_sites] spatial_deviation_t_i = rep_vector(0, n_sites);
  // To explain spatial_deviation_t_i = 0 initialization,
  // we need a spatial_deviation_t_i value for t_1 = 1, so we will
  // initialize the process with index 1 of spatial
  // deviation noise matrix w/out spatial_deviation ar coeff, so we set
  // spatial_deviation_t_i to 0 initially.

  for (t_i in 1:n_time) {
    spatial_deviation_t_i = (spatial_deviation_ar_coeff*spatial_deviation_t_i)
                            + spatial_deviation_noise_matrix[,t_i];
    log_site_rt[,t_i] = log_state_rt[t_i] + spatial_deviation_t_i;
  }

  return log_site_rt;
}
