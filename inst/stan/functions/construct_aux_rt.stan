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
vector construct_aux_rt(vector log_state_rt, real state_deviation_ar_coeff, real state_deviation_init, vector state_deviation_noise_vec) {


  // presets
  int n_time = dims(state_deviation_noise_vec)[1];
  vector[n_time] log_aux_site_rt;
  vector[n_time] state_deviation;

  state_deviation[1] = state_deviation_init;
  for (t_i in 2:n_time) {
    state_deviation[t_i] = (state_deviation_ar_coeff*state_deviation[t_i-1])
                              + state_deviation_noise_vec[t_i];
  }

  for (t_i in 1:n_time) {
    log_aux_site_rt[t_i] = log_state_rt[t_i] + state_deviation[t_i];
  }

  return log_aux_site_rt;
}
