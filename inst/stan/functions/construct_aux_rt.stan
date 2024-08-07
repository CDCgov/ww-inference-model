/**
  * Generates the aux site Rt given all the random components.
  * @param log_state_rt "State" level unadjusted Rt, on log scale.
  * @param state_deviation_ar_coeff Coefficient for AR(1) temporal correlation on
  * subpopulation deviations.
  * @param state_deviation_init Initial value of the spatial deviation
  * @param state_deviation_noise_vec vector with random vectors of noise
  * for state devaiation process for aux site.
  * @return A vector for unadjusted Rt for aux site where columns are time points.
  */
vector construct_aux_rt(vector log_state_rt,
                        real state_deviation_ar_coeff,
                        vector state_deviation_noise_vec) {


  // presets
  int n_time = dims(state_deviation_noise_vec)[1];
  vector[n_time] log_aux_site_rt = ar1(
    log_state_rt,
    state_deviation_ar_coeff,
    1,
    state_deviation_noise_vec,
    0
  );

  return log_aux_site_rt;
}
