/**
  * Generates a matrix where columns are random vectors for spatial process.
  * @param sigma_matrix Matrix of covariances, constructed externally.
  * @param n_time timepoints for number of random vectors of spatial process desired.
  * @return A matrix for spatial deviation noise in use of spatial Rt process.
  */
matrix spatial_deviation_noise_matrix_rng(matrix sigma_matrix,
                                          int n_time) {

  int n_sites = rows(sigma_matrix);
  matrix[n_sites,n_time] spatial_deviation_noise;

    for (t_i in 1:n_time) {
      spatial_deviation_noise[,t_i] = multi_normal_rng(rep_vector(0,n_sites), sigma_matrix);
    }

  return spatial_deviation_noise;
}
