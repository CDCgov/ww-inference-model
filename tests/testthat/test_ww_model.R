test_that("Test the wastewater inference model on simulated data.", {
  #######
  # run model briefly on the simulated data
  #######
  withr::with_seed(5, {
    fit <- do.call(
      silent_wwinference,
      model_test_data
    )
  })

  params <- model_test_data$model_spec$params
  obs_last_draw <- posterior::subset_draws(fit$fit$result$draws(),
    draw = 25
  )

  # Check all parameters (ignoring their dimensions) are in both fits
  # But in a way that makes error messages easy to understand
  obs_par_names <- get_nonmatrix_names_from_draws(obs_last_draw)
  exp_par_names <- get_nonmatrix_names_from_draws(test_fit_last_draw)

  expect_true(
    all(!!obs_par_names %in% !!exp_par_names)
  )

  expect_true(
    all(!!exp_par_names %in% !!obs_par_names)
  )

  # Check dims
  obs_par_lens <- get_par_dims_flat(obs_last_draw)
  exp_par_lens <- get_par_dims_flat(test_fit_last_draw)

  agg_names <- c(names(obs_par_lens), names(exp_par_lens)) |> unique()
  for (param in agg_names) {
    expect_equal(
      obs_par_lens[!!param],
      exp_par_lens[!!param]
    )
  }
  expect_mapequal(
    obs_par_lens,
    exp_par_lens
  )

  # Check the parameters we care most about
  model_params <- c(
    "eta_sd", "autoreg_rt", "log_r_mu_intercept", "sigma_rt",
    "autoreg_rt_site", "i0_over_n", "sigma_i0", "sigma_growth",
    "initial_growth", "inv_sqrt_phi_h", "sigma_ww_site_mean",
    "sigma_ww_site_sd",
    "p_hosp_w_sd", "t_peak", "dur_shed", "ww_site_mod_sd", "rt", "rt_site_t",
    "p_hosp", "w", "hosp_wday_effect", "eta_i0", "eta_growth",
    "infection_feedback", "p_hosp_mean"
  )

  for (param in model_params) {
    # Compare everything, with numerical tolerance
    testthat::expect_equal(
      obs_last_draw,
      test_fit_last_draw,
      tolerance = 0.0001
    )
  }

  # Testing draws
  model_draws <- get_draws(fit)
  expect_length(model_draws, 4)

  expect_error(get_draws(fit, what = "something else"))

  # Getting a forecast date
  forecast_date <- model_draws$predicted_counts$date
  forecast_date <- min(forecast_date) + floor(diff(range(forecast_date)) * .75)

  # Extracting the observed data for the plots
  count_data_eval <- model_draws$predicted_counts |>
    dplyr::select(observed_value, date)

  expect_true(
    inherits(
      plot(
        model_draws,
        what = "predicted_counts",
        forecast_date = forecast_date,
        n_draws_to_plot = model_test_data$fit_opts$iter_sampling,
        count_data_eval = count_data_eval,
        count_data_eval_col_name = "observed_value"
      ),
      "ggplot"
    )
  )
  expect_true(
    inherits(
      plot(
        model_draws,
        what = "predicted_ww",
        forecast_date = forecast_date,
        n_draws_to_plot = model_test_data$fit_opts$iter_sampling
      ),
      "ggplot"
    )
  )
  expect_true(
    inherits(
      plot(
        model_draws,
        what = "global_rt",
        forecast_date = forecast_date,
        n_draws_to_plot = model_test_data$fit_opts$iter_sampling
      ),
      "ggplot"
    )
  )
  expect_true(
    inherits(
      plot(
        model_draws,
        what = "subpop_rt",
        forecast_date = forecast_date,
        n_draws_to_plot = model_test_data$fit_opts$iter_sampling
      ),
      "ggplot"
    )
  )

  expect_error(plot(model_draws, what = "something else"))
})
