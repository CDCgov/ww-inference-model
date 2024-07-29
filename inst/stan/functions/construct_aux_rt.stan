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
  vector[n_time] log_aux_site_rt;
  real state_deviation_t_i = 0;
  // To explain state_deviation_t_i = 0 initialization,
  // we need a state_deviation_t_i value for t_1 = 1, so we will
  // initialize the process with index 1 of spatial
  // deviation noise matrix w/out spatial_deviation ar coeff, so we set
  // state_deviation_t_i to 0 initially.

  for (t_i in 1:n_time) {
    state_deviation_t_i = (state_deviation_ar_coeff*state_deviation_t_i)
                              + state_deviation_noise_vec[t_i];
    log_aux_site_rt[t_i] = log_state_rt[t_i] + state_deviation_t_i;
  }

  return log_aux_site_rt;
}
