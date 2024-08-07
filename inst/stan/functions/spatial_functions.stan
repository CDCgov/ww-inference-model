/**
  * File to call all desired spatial functions for use with expose_functions()
  */
functions{
  #include exponential_decay_corr_func.stan
  #include independence_corr_func.stan
  #include construct_spatial_rt.stan
  #include spatial_deviation_noise_matrix_rng.stan
  #include spatial_rt_process_rng.stan
  #include construct_aux_rt.stan
  #include state_deviation_noise_vec_aux_rng.stan
  #include aux_site_process_rng.stan
  #include matrix_normalization.stan
}
