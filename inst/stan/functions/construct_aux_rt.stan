/**
  * Generates the aux site Rt given all the random components.
  * @param log_state_rt "State" level unadjusted Rt, on log scale.
  * @param state_deviation_ar_coeff Coefficient for AR(1) temporal correlation on
  * subpopulation deviations.
  * @param scale_factor Scaling factor for aux site process.
  * @param sigma_eps Parameter for construction of covariance matrix
  * of spatial epsilon.
  * @param z Vector with random vectors of noise on standard normal.
  * @param init_bool Boolean for making initial value stationary( 1 or 0 ).
  * @return A vector for unadjusted Rt for aux site where columns are time points.
  */
vector construct_aux_rt(vector log_state_rt,
                        real state_deviation_ar_coeff,
                        real scaling_factor,
                        real sigma_eps,
                        vector z,
                        int init_bool) {


  // presets
  int n_time = dims(z)[1];
  vector[n_time] log_aux_site_rt = ar1(
    log_state_rt,
    state_deviation_ar_coeff,
    scaling_factor * sigma_eps,
    z,
    init_bool
  );

  return log_aux_site_rt;
}
