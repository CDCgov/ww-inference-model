/**
  * Generates the subpopulation Rt trajectories for a spatially correlated model
  * @param log_state_rt "State" level unadjusted Rt, on log scale.
  * @param scaling_factor Scaling factor for aux site process.
  * @param sigma_eps Parameter for construction of covariance matrix
  * of spatial epsilon.
  * @param state_deviation_ar_coeff Coefficient for AR(1) temporal correlation on
  * state deviations for aux site.
  * @param init_bool Boolean for making initial value stationary( 1 or 0 ).
  * @return A vector for unadjusted Rt for aux site where columns are time points.
  */
vector aux_site_process_rng(vector log_state_rt,
                              real scaling_factor,
                              real sigma_eps,
                              real state_deviation_ar_coeff,
                              int init_bool) {


  //Presets
  int n_time = size(log_state_rt);

  // creates vector of standard normal for contruct aux rt.
  vector[n_time] z = state_deviation_noise_vec_aux_rng(1,
                                                       1,
                                                       n_time);
  vector[n_time] log_aux_site_rt = construct_aux_rt(log_state_rt,
                                                    state_deviation_ar_coeff,
                                                    scaling_factor,
                                                    sigma_eps,
                                                    z,
                                                    init_bool);

  return log_aux_site_rt;
}
