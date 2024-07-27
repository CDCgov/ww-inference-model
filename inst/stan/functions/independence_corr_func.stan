/**
  * Assembles an independence correlation structure.
  * @param n_sites number of sites/subpopulations.
  * @return A correlation matrix for a spatially independent structure.
  */
matrix independence_corr_func(int n_sites) {

  return identity_matrix(n_sites);
}
