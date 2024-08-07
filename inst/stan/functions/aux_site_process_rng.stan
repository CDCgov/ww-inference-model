/**
  * Generates the subpopulation Rt trajectories for a spatially correlated model
  * @param log_state_rt "State" level unadjusted Rt, on log scale.
  * @param scaling_factor Scaling factor for aux site process.
  * @param sigma_eps Parameter for construction of covariance matrix
  * of spatial epsilon.
  * @param state_deviation_ar_coeff Coefficient for AR(1) temporal correlation on
  * state deviations for aux site.
  * @param init_stat Boolean. Should the initial value of the AR(1) be drawn
  * from the process's stationary stationary distribution (1) or from the
  * process's conditional error distribution (0)? Note that the process only has
  * a defined stationary distribution if `state_deviation_ar_coeff` < 1.
  * @return A vector for unadjusted Rt for aux site where columns are time points.
  */
vector aux_site_process_rng(vector log_state_rt,
                              real scaling_factor,
                              real sigma_eps,
                              real state_deviation_ar_coeff,
                              int init_stat) {


  //Presets
  int n_time = size(log_state_rt);

  // creates vector of standard normal for contruct aux rt.
  vector[n_time] z;
  for (i in 1:n_time){
    z[i] = std_normal_rng();
  }
  vector[n_time] log_aux_site_rt = construct_aux_rt(log_state_rt,
                                                    state_deviation_ar_coeff,
                                                    scaling_factor,
                                                    sigma_eps,
                                                    z,
                                                    init_stat);

  return log_aux_site_rt;
}
