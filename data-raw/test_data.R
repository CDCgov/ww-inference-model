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
  conc_col_name = "genome_copies_per_ml",
  lod_col_name = "lod"
)

hosp_data_preprocessed <- wwinference::preprocess_hosp_data(
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
generation_interval <- wwinference::generation_interval
inf_to_hosp <- wwinference::inf_to_hosp

# Assign infection feedback equal to the generation interval
infection_feedback_pmf <- generation_interval
model <- wwinference::compile_model()

model_spec <- wwinference::get_model_spec(
  forecast_date = forecast_date,
  calibration_time = calibration_time,
  forecast_horizon = forecast_horizon,
  generation_interval = generation_interval,
  inf_to_count_delay = inf_to_hosp,
  infection_feedback_pmf = infection_feedback_pmf
)

fit <- wwinference::wwinference(
  ww_data_to_fit,
  hosp_data_preprocessed,
  model_spec = model_spec,
  mcmc_options = wwinference::get_mcmc_options(
    n_chains = 1,
    iter_sampling = 25,
    iter_warmup = 25
  ),
  generate_initial_values = FALSE,
  compiled_model = model
)


# Create the toy stan data object for testing
toy_stan_data <- wwinference::get_stan_data(
  input_count_data = hosp_data_preprocessed,
  input_ww_data = ww_data_to_fit,
  forecast_date = model_spec$forecast_date,
  calibration_time = model_spec$calibration_time,
  forecast_horizon = model_spec$forecast_horizon,
  generation_interval = model_spec$generation_interval,
  inf_to_count_delay = model_spec$inf_to_count_delay,
  infection_feedback_pmf = model_spec$infection_feedback_pmf,
  params = model_spec$params,
  compute_likelihood = 1
)


# Generate the last draw of a very short run for testing
toy_stan_fit_last_draw <- posterior::subset_draws(fit$raw_fit_obj$draws(),
  draw = 25
)
# Save the data as internal data. Every time the model changes, will need
# to regenerate this testing data.
usethis::use_data(
  toy_stan_data,
  toy_stan_fit_last_draw,
  internal = TRUE,
  overwrite = TRUE
)
