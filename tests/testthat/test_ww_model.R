test_that("Test the wastewater inference model on simulated data.", {
  #######
  # run model briefly on the simulated data
  #######
  model <- compiled_site_inf_model
  fit <- model$sample(
    data = toy_stan_data,
    seed = 123,
    iter_sampling = 25,
    iter_warmup = 25,
    chains = 1
  )

  obs_last_draw <- posterior::subset_draws(fit$draws(), draw = 25)

  # Check all parameters (ignoring their dimensions) are in both fits
  # But in a way that makes error messages easy to understand
  obs_par_names <- get_nonmatrix_names_from_draws(obs_last_draw)
  exp_par_names <- get_nonmatrix_names_from_draws(toy_stan_fit_last_draw)

  expect_true(
    all(!!obs_par_names %in% !!exp_par_names)
  )

  expect_true(
    all(!!exp_par_names %in% !!obs_par_names)
  )

  # Check dims
  obs_par_lens <- get_par_dims_flat(obs_last_draw)
  exp_par_lens <- get_par_dims_flat(toy_stan_fit_last_draw)

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
      toy_stan_fit_last_draw,
      tolerance = 0.0001
    )
  }
})
