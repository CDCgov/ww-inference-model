############
# Make entirely fake stan input data via prior-predictive generated quantities
############

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
  iter_warmup = 25,
  iter_sampling = 25,
  chains = 1,
  show_messages = FALSE
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

withr::with_seed(5, {
  fit <- do.call(
    wwinference::wwinference,
    model_test_data
  )
})


# Generate the last draw of a very short run for testing
test_fit_last_draw <- posterior::subset_draws(
  fit$fit$result$draws(),
  draw = 25
)
# Save the data as internal data. Every time the model changes, will need
# to regenerate this testing data.
usethis::use_data(
  model_test_data,
  test_fit_last_draw,
  internal = TRUE,
  overwrite = TRUE
)
