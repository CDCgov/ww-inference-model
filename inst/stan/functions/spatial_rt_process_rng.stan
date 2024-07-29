/**
  * Generates the subpopulation Rt trajectories for a spatially correlated model
  * @param log_state_rt "State" level unadjusted Rt, on log scale.
  * @param Sigma_matrix Matrix of covariances, constructed externally.
  * @param spatial_deviation_ar_coeff Coefficient for AR(1) temporal correlation on
  * subpopulation deviations.
  * @return A matrix for unadjusted Rt where rows are subpopulations
  * and columns are time points.
  */
matrix spatial_rt_process_rng(vector log_state_rt,
                          matrix Sigma_matrix,
                          real spatial_deviation_ar_coeff) {


  //Presets
  int n_time = dims(log_state_rt)[1];
  int n_sites = rows(Sigma_matrix);

  matrix[n_sites,n_time] spatial_deviation_noise_matrix = spatial_deviation_noise_matrix_rng(Sigma_matrix,
                                                                                             n_time);
  matrix[n_sites,n_time] log_site_rt = construct_spatial_rt(log_state_rt,
                                                            spatial_deviation_ar_coeff,
                                                            spatial_deviation_noise_matrix);

  return log_site_rt;
}
