/**
  * Generates the subpopulation Rt trajectories for a spatially correlated model
  * given all the random components.
  * @param log_state_rt "State" level unadjusted Rt, on log scale.
  * @param spatial_deviation_ar_coeff Coefficient for AR(1) temporal correlation on
  * subpopulation deviations.
  * @param spatial_deviation_init Initial vector of the spatial deviation
  * @param spatial_deviation_noise_matrix matrix with random vectors of noise
  * for spatial devaiation process.
  * @return A matrix for unadjusted Rt where rows are subpopulations
  * and columns are time points.
  */
matrix construct_spatial_rt(vector log_state_rt,
                            real spatial_deviation_ar_coeff,
                            vector spatial_deviation_init,
                            matrix spatial_deviation_noise_matrix) {


  // presets
  int n_time = cols(spatial_deviation_noise_matrix);
  int n_sites = rows(spatial_deviation_noise_matrix);
  matrix[n_sites,n_time] log_site_rt;
  matrix[n_sites,n_time] spatial_deviation;

  spatial_deviation[,1] = spatial_deviation_init;
  for (t_i in 2:n_time) {
    spatial_deviation[,t_i] = (spatial_deviation_ar_coeff*spatial_deviation[,t_i-1])
                              + spatial_deviation_noise_matrix[,t_i];
  }

  for (t_i in 1:n_time) {
    log_site_rt[,t_i] = (log_state_rt[t_i]*rep_vector(1,n_sites))
                        + spatial_deviation[,t_i];
  }

  return log_site_rt;
}
