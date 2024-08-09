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
inf_to_count_delay <- to_simplex(c(rep(0.01, 12), rep(0.2, 4), rep(0.01, 10)))
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


test_that("Test that function returns list with the necessary elements", {
  result <- stan_data_list <- get_stan_data(
    input_count_data,
    input_ww_data,
    forecast_date,
    forecast_horizon,
    calibration_time,
    generation_interval,
    inf_to_count_delay,
    infection_feedback_pmf,
    params,
    include_ww
  )


  testthat::expect_true(all(names(result) %in% c("input_data", "stan_args")))
})
