/**
  * Assembles the correlation matrix with an exponential decay corr. func.
  * @param dist_matrix matrix with distances.
  * @param phi scalar parameter controlling the decay rate of correlation with
  * distance. Should be non-negative. Smaller values imply a greater reduction
  * in correlation with increased distance.
  * @param l parameter for function.
  * @return A matrix, the correlation matrix for exp. corr. func.
  */
matrix exponential_decay_corr_func( matrix dist_matrix, real phi, real l) {

  int k = rows(dist_matrix);
  matrix[k, k] corr_matrx = exp(-pow(dist_matrix / phi, l));

  return corr_matrx;
}
