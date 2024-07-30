/**
  * Generates a vector for aux site's state deviation noise.
  * @param scaling_factor Scaling factor for aux site process.
  * @param sigma_eps Parameter for construction of covariance matrix
  * of spatial epsilon.
  * @param n_time timepoints for number of random vectors of spatial process desired.
  * @return A vector for state deviation noise in use of auz site Rt process.
  */
vector state_deviation_noise_vec_aux_rng(real scaling_factor,
                                         real sigma_eps,
                                         int n_time) {

  vector[n_time] state_deviation_noise = to_vector(
    normal_rng(rep_vector(0, n_time), scaling_factor * sigma_eps)
    );

  return state_deviation_noise;
}
