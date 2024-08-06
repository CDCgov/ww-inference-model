test_that(
  "Test if stan functions for AUX site is correct",
  {
    model <- compiled_site_inf_model
    space_model_fxns <- spatial_fxns$functions

    log_state_rt <- rnorm(
      n = 1000,
      mean = 1.2,
      sd = 0.05
    )
    scaling_factor <- 1.1
    sigma_eps <- sqrt(0.02)
    state_deviation_ar_coeff <- 0.1
    n_times <- length(log_state_rt)


    passed_tests <- 0
    num_tests <- 100

    state_deviation_noise_vec <- rnorm(
      n = n_times,
      mean = 0,
      sd = scaling_factor * sigma_eps
    )
    stan_log_site_rt <- space_model_fxns$construct_aux_rt(
      log_state_rt,
      state_deviation_ar_coeff,
      state_deviation_noise_vec
    )

    r_log_aux_site_rt <- matrix(
      data = 0,
      nrow = 1,
      ncol = n_times
    )
    state_deviation_t_i <- 0
    for (t_i in 1:n_times) {
      state_deviation_t_i <- state_deviation_ar_coeff * state_deviation_t_i +
        state_deviation_noise_vec[t_i]
      r_log_aux_site_rt[t_i] <- log_state_rt[t_i] + state_deviation_t_i
    }
    r_log_aux_site_rt <- as.vector(r_log_aux_site_rt)


    testthat::expect_equal(
      stan_log_site_rt,
      r_log_aux_site_rt,
      ignore_attr = TRUE
    )
  }
)
