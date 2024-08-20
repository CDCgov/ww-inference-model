/**
  * Normalizes a distance matrix using its largest distance.
  * @param dist_matrx A distance matrix.
  * @return A distance matrix that has been normalized.
  *
  */
matrix dist_matrix_normalization(matrix dist_matrx) {
  int n = cols(dist_matrx);
  real max_val = max(dist_matrx);
  matrix[n,n] norm_dist_matrx = dist_matrx / max_val;
  return norm_dist_matrx;
}
