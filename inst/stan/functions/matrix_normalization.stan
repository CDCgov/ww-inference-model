/**
  * Normalizes a matrix using its determinant to the power of 1 / n.
  * @param matrx A matrix, usually the correlation matrix need to normalize
  * @return A matrix that has been normalized.
  *
  */
matrix matrix_normalization(matrix matrx) {
  int n = cols(matrx);
  matrix[n,n] norm_matrx = matrx / pow(determinant(matrx), 1.0 / n);
  return norm_matrx;
}
