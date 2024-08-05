test_that(
  "Test if non-centered parameter transformation of MVN(0,I) agrees with stan
  spatial_deviation_noise_matrix_rng.stan",
  {
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

    passed_tests <- 0
    num_tests <- 100
    for (i in 1:num_tests) {
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
      chol_sigma_matrix <- t(chol(
        sigma_matrix
      ))
      r_noise_matrix <- apply(
        t(r_non_center_noise_matrix),
        2,
        function(x) chol_sigma_matrix %*% x
      )

      # comparing with cramer test
      cramer_p_value <- cramer.test(
        t(stan_noise_matrix),
        t(r_noise_matrix)
      )$p.value

      # updating passed tests
      if (cramer_p_value > 0.01) {
        passed_tests <- passed_tests + 1
      }
    }

    testthat::expect_gt(
      passed_tests,
      num_tests * .95
    )
  }
)
