# Generate test data
ww_data <- tibble::tibble(
  date = rep(seq(
    from = lubridate::ymd("2023-08-01"),
    to = lubridate::ymd("2023-11-01"),
    by = "weeks"
  ), 2),
  site = c(rep(1, 14), rep(2, 14)),
  lab = c(rep(1, 28)),
  conc = abs(rnorm(28, mean = 500, sd = 50)),
  lod = c(rep(20, 14), rep(15, 14)),
  site_pop = c(rep(2e5, 14), rep(4e5, 14))
)

ww_data_preprocessed <- preprocess_ww_data(ww_data,
  conc_col_name = "conc",
  lod_col_name = "lod"
)
input_ww_data <- indicate_ww_exclusions(ww_data_preprocessed)

hosp_data <- tibble::tibble(
  date = seq(
    from = lubridate::ymd("2023-07-01"),
    to = lubridate::ymd("2023-10-30"),
    by = "days"
  ),
  daily_admits = sample(5:70, 122, replace = TRUE),
  state_pop = rep(1e6, 122)
)

input_count_data <- preprocess_count_data(
  hosp_data,
  "daily_admits",
  "state_pop"
)

generation_interval <- to_simplex(c(0.01, 0.2, 0.3, 0.2, 0.1, 0.1, 0.01))
inf_to_count_delay <- to_simplex(c(
  rep(0.01, 12), rep(0.2, 4),
  rep(0.01, 10)
))
infection_feedback_pmf <- generation_interval

params <- get_params(
  system.file("extdata", "example_params.toml",
    package = "wwinference"
  )
)
forecast_date <- "2023-11-06"
calibration_time <- 90
forecast_horizon <- 28
include_ww <- 1


test_that("wwinference model can compile", {
  expect_no_error(compile_model())
})

test_that("Function to get mcmc options produces the expected outputs", {
  mcmc_options <- get_mcmc_options()
  expected_names <- c(
    "iter_warmup", "iter_sampling", "seed", "adapt_delta", "max_treedepth"
  )
  checkmate::expect_names(names(mcmc_options), must.include = expected_names)
})

test_that("Function to get model specs produces expected outputs", {
  model_spec <- get_model_spec()
  expected_names <- c(
    "generation_interval", "inf_to_count_delay",
    "infection_feedback_pmf", "include_ww",
    "compute_likelihood", "params"
  )
  # Checkmade doesn't work here for a list, says it must be a character vector
  expect_true(all(names(model_spec) %in% expected_names))
})

test_that("Passing invalid args to fit_opts throws an error ", {
  expect_error(
    wwinference(
      ww_data = input_ww_data,
      count_data = input_count_data,
      forecast_date = forecast_date,
      model_spec = get_model_spec,
      fit_opts = list(not_an_arg = 4)
    ),
    regexp = c("Names must be a subset of ")
  )
})
