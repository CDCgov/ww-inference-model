/**
  * Generates the subpopulation Rt trajectories for a spatially correlated model
  * @param log_state_rt "State" level unadjusted Rt, on log scale.
  * @param dist_matrix Distance matrix for use in correaltion function.
  * @param sigma_eps Parameter for construction of covariance matrix
  * of spatial epsilon.
  * @param phi Parameter for construction of correaltion matrix.
  * @param l Parameter for construction of correlation matrix.
  * @param spatial_deviation_ar_coeff Coefficient for AR(1) temporal correlation on
  * subpopulation deviations.
  * @param spatial_deviation_init Initial vector of spatial deviation.
  * @return A matrix for unadjusted Rt where rows are subpopulations
  * and columns are time points.
  */
matrix spatial_rt_process_rng(vector log_state_rt,
                          matrix dist_matrix,
                          real sigma_eps,
                          real phi,
                          real l,
                          real spatial_deviation_ar_coeff,
                          vector spatial_deviation_init) {


  //Presets
  int n_time = dims(log_state_rt)[1];
  int n_sites = rows(dist_matrix);

  matrix[n_sites,n_time] spatial_deviation_noise_matrix = spatial_deviation_noise_matrix_rng(dist_matrix,
                                                                                             sigma_eps,
                                                                                             phi,
                                                                                             l,
                                                                                             n_time);
  matrix[n_sites,n_time] log_site_rt = construct_spatial_rt(log_state_rt,
                                                            spatial_deviation_ar_coeff,
                                                            spatial_deviation_init,
                                                            spatial_deviation_noise_matrix);

  return log_site_rt;
}
