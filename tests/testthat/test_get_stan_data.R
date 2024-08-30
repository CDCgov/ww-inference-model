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
ww_data_filtered <- indicate_ww_exclusions(ww_data_preprocessed)


hosp_data <- tibble::tibble(
  date = seq(
    from = lubridate::ymd("2023-07-01"),
    to = lubridate::ymd("2023-10-30"),
    by = "days"
  ),
  daily_admits = sample(5:70, 122, replace = TRUE),
  state_pop = rep(1e6, 122)
)

count_data <- preprocess_count_data(
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

input_count_data <- get_input_count_data_for_stan(
  count_data,
  calibration_time
)
first_count_data_date <- min(input_count_data$date, na.rm = TRUE)
last_count_data_date <- max(input_count_data$date, na.rm = TRUE)
input_ww_data <- get_input_ww_data_for_stan(
  ww_data_filtered,
  first_count_data_date,
  last_count_data_date,
  calibration_time
)


test_that(paste0(
  "Test that modifying calibration time generates data of expected",
  " length"
), {
  result <- get_input_count_data_for_stan(
    count_data,
    calibration_time = 80
  )
  expect_true(nrow(result) == 80)
})

test_that(paste0(
  "Test that things not flagged for removal don't get removed ",
  "and things that are flagged for removal do get removed"
), {
  ww_data_no_exclusions <- ww_data_filtered
  ww_data_no_exclusions$exclude <- 0
  input_ww_data_ne <- get_input_ww_data_for_stan(
    ww_data_no_exclusions,
    first_count_data_date,
    last_count_data_date,
    calibration_time
  )

  expect_true(nrow(input_ww_data_ne) == nrow(input_ww_data))

  ww_data_w_exclusions <- ww_data_filtered
  ww_data_w_exclusions$exclude[10] <- 1
  input_ww_data_we <- get_input_ww_data_for_stan(
    ww_data_w_exclusions,
    first_count_data_date,
    last_count_data_date,
    calibration_time
  )

  expect_true(nrow(input_ww_data_ne) == nrow(input_ww_data_we) + 1)
})



test_that(paste0(
  "Test that passing out of window wastewater data behaves as",
  "expected"
), {
  # Make wastewater data outside of scope of admissions data
  recent_ww_data <- tibble::tibble(
    date = rep(seq(
      from = lubridate::ymd("2024-08-01"),
      to = lubridate::ymd("2024-11-01"),
      by = "weeks"
    ), 2),
    site = c(rep(1, 14), rep(2, 14)),
    lab = c(rep(1, 28)),
    conc = abs(rnorm(28, mean = 500, sd = 50)),
    lod = c(rep(20, 14), rep(15, 14)),
    site_pop = c(rep(2e5, 14), rep(4e5, 14))
  )

  recent_ww_data_preprocessed <- preprocess_ww_data(recent_ww_data,
    conc_col_name = "conc",
    lod_col_name = "lod"
  )
  recent_input_ww_data <- indicate_ww_exclusions(recent_ww_data_preprocessed)

  recent_input_ww_data_for_stan <- get_input_ww_data_for_stan(
    recent_input_ww_data,
    first_count_data_date,
    last_count_data_date,
    calibration_time
  )

  expect_error(get_stan_data(
    input_count_data,
    recent_input_ww_data_for_stan,
    forecast_date,
    forecast_horizon,
    calibration_time,
    generation_interval,
    inf_to_count_delay,
    infection_feedback_pmf,
    params,
    include_ww,
    dist_matrix = NULL,
    corr_structure_switch = 0
  ))

  # Make wastewater data outside of scope of admissions data
  old_ww_data <- tibble::tibble(
    date = rep(seq(
      from = lubridate::ymd("2022-08-01"),
      to = lubridate::ymd("2022-11-01"),
      by = "weeks"
    ), 2),
    site = c(rep(1, 14), rep(2, 14)),
    lab = c(rep(1, 28)),
    conc = abs(rnorm(28, mean = 500, sd = 50)),
    lod = c(rep(20, 14), rep(15, 14)),
    site_pop = c(rep(2e5, 14), rep(4e5, 14))
  )

  old_ww_data_preprocessed <- preprocess_ww_data(old_ww_data,
    conc_col_name = "conc",
    lod_col_name = "lod"
  )
  old_input_ww_data <- indicate_ww_exclusions(old_ww_data_preprocessed)

  old_input_ww_data_for_stan <- get_input_ww_data_for_stan(
    old_input_ww_data,
    first_count_data_date,
    last_count_data_date,
    calibration_time
  )


  expect_error(get_stan_data(
    input_count_data,
    old_input_ww_data,
    forecast_date,
    forecast_horizon,
    calibration_time,
    generation_interval,
    inf_to_count_delay,
    infection_feedback_pmf,
    params,
    include_ww,
    dist_matrix = NULL,
    corr_structure_switch = 0
  ))
})

test_that("Test that pmf check works as expected", {
  expect_warning(get_stan_data(
    input_count_data,
    input_ww_data,
    forecast_date,
    forecast_horizon,
    calibration_time,
    generation_interval = to_simplex(rep(1, 100)),
    inf_to_count_delay,
    infection_feedback_pmf,
    params,
    include_ww,
    dist_matrix = NULL,
    corr_structure_switch = 0
  ))

  expect_warning(get_stan_data(
    input_count_data,
    input_ww_data,
    forecast_date,
    forecast_horizon,
    calibration_time,
    generation_interval,
    inf_to_count_delay = to_simplex(rep(1, 100)),
    infection_feedback_pmf,
    params,
    include_ww,
    dist_matrix = NULL,
    corr_structure_switch = 0
  ))

  expect_warning(get_stan_data(
    input_count_data,
    input_ww_data,
    forecast_date,
    forecast_horizon,
    calibration_time,
    generation_interval,
    inf_to_count_delay,
    infection_feedback_pmf = to_simplex(rep(1, 100)),
    params,
    include_ww,
    dist_matrix = NULL,
    corr_structure_switch = 0
  ))

  expect_error(get_stan_data(
    input_count_data,
    input_ww_data,
    forecast_date,
    forecast_horizon,
    calibration_time,
    generation_interval,
    inf_to_count_delay,
    infection_feedback_pmf = c(0.5, 0.4, 0.2),
    params,
    include_ww,
    dist_matrix = NULL,
    corr_structure_switch = 0
  ))
})
