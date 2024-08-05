/**
  * Generates the subpopulation Rt trajectories for a spatially correlated model
  * @param log_state_rt "State" level unadjusted Rt, on log scale.
  * @param scaling_factor Scaling factor for aux site process.
  * @param sigma_eps Parameter for construction of covariance matrix
  * of spatial epsilon.
  * @param state_deviation_ar_coeff Coefficient for AR(1) temporal correlation on
  * state deviations for aux site.
  * @return A vector for unadjusted Rt for aux site where columns are time points.
  */
vector aux_site_process_rng(vector log_state_rt,
                              real scaling_factor,
                              real sigma_eps,
                              real state_deviation_ar_coeff) {


  //Presets
  int n_time = size(log_state_rt);

  vector[n_time] state_deviation_noise_vec = state_deviation_noise_vec_aux_rng(scaling_factor,
                                                                               sigma_eps,
                                                                               n_time);
  vector[n_time] log_aux_site_rt = construct_aux_rt(log_state_rt,
                                                    state_deviation_ar_coeff,
                                                    state_deviation_noise_vec);

  return log_aux_site_rt;
}
