/**
  * Generates a matrix where columns are random vectors for spatial process.
  * Currently this is for an exponential correlation function.
  * @param dist_matrix Distance matrix for correlation constr.
  * @param sigma_eps Parameter for construction of covariance matrix
  * of spatial epsilon.
  * @param phi Phi parameter for exponential decay corr. func.
  * @param l l paraamter for exponential decay corr. func.
  * @param n_time timepoints for number of random vectors of spatial process desired.
  * @return A matrix for spatial deviation noise in use of spatial Rt process.
  */
matrix spatial_deviation_noise_matrix_rng(matrix dist_matrix,
                                          real sigma_eps,
                                          real phi,
                                          real l,
                                          int n_time) {

  int n_sites = rows(dist_matrix);
  // correlation matrix constr.
  matrix[n_sites,n_sites] omega_matrix = exponential_decay_corr_func(dist_matrix, phi, l);
  matrix[n_sites,n_sites] sigma_matrix = sigma_eps^2 * omega_matrix;


  matrix[n_sites,n_time] spatial_deviation_noise;

    for (t_i in 1:n_time) {
      spatial_deviation_noise[,t_i] = multi_normal_rng(rep_vector(0,n_sites), sigma_matrix);
    }

  return spatial_deviation_noise;
}
