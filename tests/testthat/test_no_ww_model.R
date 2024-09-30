options(cmdstanr_warn_inits = FALSE)

hosp_data <- wwinference::hosp_data
ww_data <- wwinference::ww_data
params <- wwinference::get_params(
  fs::path_package("extdata", "example_params.toml",
    package = "wwinference"
  )
)


# Data pre-processing --------------------------------------------------------
ww_data_preprocessed <- wwinference::preprocess_ww_data(
  ww_data,
  conc_col_name = "log_genome_copies_per_ml",
  lod_col_name = "log_lod"
)

hosp_data_preprocessed <- wwinference::preprocess_count_data(
  hosp_data,
  count_col_name = "daily_hosp_admits",
  pop_size_col_name = "state_pop"
)

ww_data_to_fit <- wwinference::indicate_ww_exclusions(
  ww_data_preprocessed,
  outlier_col_name = "flag_as_ww_outlier",
  remove_outliers = TRUE
)

forecast_date <- "2023-12-06"
calibration_time <- 90
forecast_horizon <- 28
generation_interval <- wwinference::default_covid_gi
inf_to_hosp <- wwinference::default_covid_inf_to_hosp

# Assign infection feedback equal to the generation interval
infection_feedback_pmf <- generation_interval

model_spec <- wwinference::get_model_spec(
  generation_interval = generation_interval,
  inf_to_count_delay = inf_to_hosp,
  infection_feedback_pmf = infection_feedback_pmf,
  params = params
)

mcmc_options <- list(
  seed = 5,
  iter_warmup = 500,
  iter_sampling = 250,
  chains = 2,
  show_messages = FALSE,
  show_exceptions = FALSE,
)

generate_initial_values <- TRUE

model_test_data <- list(
  ww_data = ww_data_to_fit,
  count_data = hosp_data_preprocessed,
  forecast_date = forecast_date,
  calibration_time = calibration_time,
  forecast_horizon = forecast_horizon,
  model_spec = model_spec,
  fit_opts = mcmc_options,
  generate_initial_values = generate_initial_values
)


test_that("Test that the model runs on simulated data when include_ww=0.", {
  #######
  # run model briefly on the simulated data
  #######
  model_test_data_no_ww <- model_test_data
  model_test_data_no_ww$model_spec$include_ww <- 0
  model_test_data_no_ww$compiled_model <- compiled_site_inf_model

  expect_no_error(withr::with_seed(55, {
    fit <- do.call(
      wwinference::wwinference,
      model_test_data_no_ww
    )
  }))
})

test_that("Test that the model runs without wastewater, include_ww=0.", {
  #######
  # run model briefly on the simulated data
  #######
  model_test_data_no_ww <- model_test_data
  model_test_data_no_ww$model_spec$include_ww <- 0
  model_test_data_no_ww$ww_data <- tibble::tibble()
  model_test_data_no_ww$compiled_model <- compiled_site_inf_model

  expect_no_error(withr::with_seed(55, {
    fit <- do.call(
      wwinference::wwinference,
      model_test_data_no_ww
    )
  }))
})

test_that("Test that the model runs without wastewater, include_ww=1.", {
  #######
  # run model briefly on the simulated data
  #######
  model_test_data_no_ww <- model_test_data
  model_test_data_no_ww$model_spec$include_ww <- 1
  model_test_data_no_ww$ww_data <- tibble::tibble()
  model_test_data_no_ww$compiled_model <- compiled_site_inf_model

  expect_no_error(withr::with_seed(55, {
    fit <- do.call(
      wwinference::wwinference,
      model_test_data_no_ww
    )
  }))
})
