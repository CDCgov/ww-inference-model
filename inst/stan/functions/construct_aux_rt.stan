/**
  * Generates the aux site Rt given all the random components.
  * @param log_state_rt "State" level unadjusted Rt, on log scale.
  * @param state_deviation_ar_coeff Coefficient for AR(1) temporal correlation on
  * subpopulation deviations.
  * @param scale_factor Scaling factor for aux site process.
  * @param init_stat Boolean. Should the initial value of the AR(1) be drawn
  * from the process's stationary stationary distribution (1) or from the
  * process's conditional error distribution (0)? Note that the process only has
  * a defined stationary distribution if `state_deviation_ar_coeff` < 1.
  * @return A vector for unadjusted Rt for aux site where columns are time points.
  */
vector construct_aux_rt(vector log_state_rt,
                        real state_deviation_ar_coeff,
                        real scaling_factor,
                        real sigma_eps,
                        vector z,
                        int init_stat) {


  // presets
  int n_time = dims(z)[1];
  vector[n_time] log_aux_site_rt = ar1(
    log_state_rt,
    state_deviation_ar_coeff,
    scaling_factor * sigma_eps,
    z,
    init_stat
  );

  return log_aux_site_rt;
}
