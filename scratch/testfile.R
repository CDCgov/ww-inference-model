# Expose the stan functions into the global environment
model <- cmdstanr::cmdstan_model(
  stan_file = file.path("inst", "stan", "functions", "spatial_functions.stan"),
  compile = TRUE,
  compile_standalone = TRUE,
  force_recompile = TRUE
)
model$expose_functions(global = TRUE)


dist_matrix1 <- as.matrix(
  dist(
    data.frame(
      x = c(85, 37, 48, 7),
      y = c(12, 75, 81, 96),
      diag = TRUE,
      upper = TRUE
    )
  )
)
dist_matrix2 <- as.matrix(
  dist(
    data.frame(
      x = c(850, 370, 480, 70),
      y = c(120, 750, 810, 960),
      diag = TRUE,
      upper = TRUE
    )
  )
)

corr_matrix1 <- exponential_decay_corr_func(
  dist_matrix = dist_matrix,
  phi = 25,
  l = 1
)
corr_matrix2 <- exponential_decay_corr_func(
  dist_matrix = dist_matrix2,
  phi = 250,
  l = 1
)
matrix_normalization(
  corr_matrix
)
