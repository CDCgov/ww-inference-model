/**
  * Generates a vector for state deviation noise for aux site.
  * Currently this is for an exponential correaltion function.
  * @param scaling_factor Scaling factor for aux site process.
  * @param sigma_eps Parameter for construction of covariance matrix
  * of spatial epsilon.
  * @param n_time timepoints for number of random vectors of spatial process desired.
  * @return A vector for state deviation noise in use of auz site Rt process.
  */
vector state_deviation_noise_vec_aux_rng(real scaling_factor,
                                         real sigma_eps,
                                         int n_time) {


  vector[n_time] state_deviation_noise;

    for (t_i in 1:n_time) {
      state_deviation_noise[t_i] = normal_rng(0, scaling_factor * sigma_eps);

  return state_deviation_noise;
}
