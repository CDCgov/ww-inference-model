test_that(
  "Test Spatial Func in Stan agrees with R",
  {
    set.seed(2024)
    model <- compiled_site_inf_model
    space_model_fxns <- spatial_fxns

    log_state_rt <- rnorm(
      n = 1000,
      mean = 1.2,
      sd = 0.02
    )
    corr_function <- exponential_decay_corr_func_r
    corr_func_params <- list(
      dist_matrix = as.matrix(dist(
        data.frame(
          x = runif(5, 0, 100),
          y = runif(5, 0, 100)
        ),
        diag = TRUE,
        upper = TRUE
      )),
      phi = 25,
      l = 1
    )
    phi_rt <- 0.8
    sigma_eps <- sqrt(0.02)
    sigma_matrix <- sigma_eps^2 * corr_function(corr_func_params)

    log_site_rt_stan <- space_model_fxns$functions$spatial_rt_process_rng(
      log_state_rt = log_state_rt,
      sigma_matrix = sigma_matrix,
      spatial_deviation_ar_coeff = phi_rt
    )
    log_site_rt_r <- spatial_rt_process(
      log_state_rt = log_state_rt,
      corr_function = corr_function,
      corr_function_params = corr_func_params,
      phi_rt = phi_rt,
      sigma_eps = sigma_eps
    )
    # these functions creates the spatial Rt for the site level given
    # the correlation or covariance structure.
    # We use cramer test for multi-variate random vectors

    cramer_p_value <- cramer.test(
      log_site_rt_r,
      log_site_rt_stan
    )$p.value

    testthat::expect_gt(
      cramer_p_value,
      0.05
    )
  }
)
