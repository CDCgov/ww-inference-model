test_that(
  "Test if non-centered parameter transformation of MVN(0,I) agrees with stan
  spatial_deviation_noise_matrix_rng.stan",
  {
    set.seed(2024)
    model <- compiled_site_inf_model
    space_model_fxns <- spatial_fxns$functions

    corr_fun_params <- list(
      dist_matrix = as.matrix(
        dist(
          data.frame(
            x = c(85, 37, 48, 7),
            y = c(12, 75, 81, 96),
            diag = TRUE,
            upper = TRUE
          )
        )
      ),
      phi = 25,
      l = 1
    )
    omega_matrix <- space_model_fxns$exponential_decay_corr_func(
      corr_fun_params$dist_matrix,
      corr_fun_params$phi,
      corr_fun_params$l
    )
    sigma_matrix <- sqrt(0.02) * omega_matrix
    n_time <- 150

    stan_noise_matrix <- space_model_fxns$spatial_deviation_noise_matrix_rng(
      sigma_matrix,
      n_time
    )


    # non-centered param
    n_sites <- ncol(
      corr_fun_params$dist_matrix
    )
    r_non_center_noise_matrix <- mvrnorm(
      n = n_time,
      mu = rep(0, n_sites),
      Sigma = diag(n_sites)
    )
    sqrt_sigma_matrix <- expm::sqrtm(
      sigma_matrix
    )
    r_noise_matrix <- apply(
      t(r_non_center_noise_matrix),
      2,
      function(x) sqrt_sigma_matrix %*% x
    )

    # comparing with cramer test
    cramer_p_value <- cramer.test(
      stan_noise_matrix,
      r_noise_matrix
    )$p.value


    testthat::expect_gt(
      cramer_p_value,
      0.05
    )
  }
)
