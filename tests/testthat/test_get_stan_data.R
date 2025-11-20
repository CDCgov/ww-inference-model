withr::with_seed(123, {
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
})


ww_data_preprocessed <- preprocess_ww_data(ww_data,
  conc_col_name = "conc",
  lod_col_name = "lod"
)
ww_data_filtered <- indicate_ww_exclusions(ww_data_preprocessed)

withr::with_seed(123, {
  hosp_data <- tibble::tibble(
    date = seq(
      from = lubridate::ymd("2023-07-01"),
      to = lubridate::ymd("2023-10-30"),
      by = "days"
    ),
    daily_admits = sample(5:70, 122, replace = TRUE),
    state_pop = rep(1e6, 122)
  )
})

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
date_time_spine <- get_date_time_spine(
  forecast_date = forecast_date,
  input_count_data = input_count_data,
  last_count_data_date = last_count_data_date,
  forecast_horizon = forecast_horizon,
  calibration_time = calibration_time
)

lab_site_site_spine <- get_lab_site_site_spine(
  input_ww_data = input_ww_data
)

site_subpop_spine <- get_site_subpop_spine(
  input_ww_data = input_ww_data,
  input_count_data = input_count_data
)

lab_site_subpop_spine <- get_lab_site_subpop_spine(
  lab_site_site_spine = lab_site_site_spine,
  site_subpop_spine = site_subpop_spine
)


test_that(paste0(
  "Test that the number of subpopulations is correct for the",
  "standard case where sum(site_pops) < total_pop"
), {
  stan_data <- get_stan_data(
    input_count_data,
    input_ww_data,
    date_time_spine,
    lab_site_site_spine,
    site_subpop_spine,
    lab_site_subpop_spine,
    last_count_data_date,
    first_count_data_date,
    forecast_date,
    forecast_horizon,
    calibration_time,
    generation_interval,
    inf_to_count_delay,
    infection_feedback_pmf,
    params,
    include_ww
  )

  expect_equal(stan_data$n_subpop, (stan_data$n_ww_sites + 1))
  expect_equal(length(stan_data$subpop_size), stan_data$n_subpops)
})

test_that(paste0(
  "Test that the number of subpopulations is correct for the ",
  "standard case where sum(site_pops) > total_pop"
), {
  input_count_data_mod <- input_count_data
  input_count_data_mod$total_pop <- sum(unique(input_ww_data$site_pop) - 100)
  site_subpop_spine_mod <- get_site_subpop_spine(
    input_ww_data = input_ww_data,
    input_count_data = input_count_data_mod
  )

  lab_site_subpop_spine_mod <- get_lab_site_subpop_spine(
    lab_site_site_spine = lab_site_site_spine,
    site_subpop_spine = site_subpop_spine_mod
  )

  expect_warning({
    stan_data_mod <- get_stan_data(
      input_count_data_mod,
      input_ww_data,
      date_time_spine,
      lab_site_site_spine,
      site_subpop_spine_mod,
      lab_site_subpop_spine_mod,
      last_count_data_date,
      first_count_data_date,
      forecast_date,
      forecast_horizon,
      calibration_time,
      generation_interval,
      inf_to_count_delay,
      infection_feedback_pmf,
      params,
      include_ww
    )
  })

  expect_equal(stan_data_mod$n_subpop, (stan_data_mod$n_ww_sites))
  expect_equal(length(stan_data_mod$subpop_size), stan_data_mod$n_subpops)
  expect_equal(stan_data_mod$norm_pop, sum(stan_data_mod$subpop_size))
})

test_that(paste0(
  "Test that the model handles include_ww = 0 ",
  "appropriately by only estimating one subpopulation"
), {
  # This happens upstream in wwinference
  input_ww_data_mod <- NULL
  site_subpop_spine_mod <- get_site_subpop_spine(
    input_ww_data = input_ww_data_mod,
    input_count_data = input_count_data
  )

  lab_site_subpop_spine_mod <- get_lab_site_subpop_spine(
    lab_site_site_spine = lab_site_site_spine,
    site_subpop_spine = site_subpop_spine_mod
  )

  stan_data_ho <- get_stan_data(
    input_count_data,
    input_ww_data_mod,
    date_time_spine,
    lab_site_site_spine,
    site_subpop_spine_mod,
    lab_site_subpop_spine_mod,
    last_count_data_date,
    first_count_data_date,
    forecast_date,
    forecast_horizon,
    calibration_time,
    generation_interval,
    inf_to_count_delay,
    infection_feedback_pmf,
    params,
    include_ww = 0
  )

  expect_equal(stan_data_ho$n_subpops, 1)
  expect_equal(length(stan_data_ho$subpop_size), 1)
})

test_that(paste0(
  "Test that the model handles include_ww = 0 ",
  "and no data appropriately"
), {
  null_ww_data <- NULL

  site_subpop_spine_mod <- get_site_subpop_spine(
    input_ww_data = null_ww_data,
    input_count_data = input_count_data
  )

  lab_site_subpop_spine_mod <- get_lab_site_subpop_spine(
    lab_site_site_spine = lab_site_site_spine,
    site_subpop_spine = site_subpop_spine_mod
  )

  stan_data_ho <- get_stan_data(
    input_count_data,
    input_ww_data = null_ww_data,
    date_time_spine,
    lab_site_site_spine,
    site_subpop_spine_mod,
    lab_site_subpop_spine_mod,
    last_count_data_date,
    first_count_data_date,
    forecast_date,
    forecast_horizon,
    calibration_time,
    generation_interval,
    inf_to_count_delay,
    infection_feedback_pmf,
    params,
    include_ww = 0
  )

  expect_equal(stan_data_ho$n_subpops, 1)
  expect_equal(length(stan_data_ho$subpop_size), 1)
})


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
  withr::with_seed(123, {
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
  })

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
  date_time_spine <- get_date_time_spine(
    forecast_date = forecast_date,
    input_count_data = input_count_data,
    last_count_data_date = last_count_data_date,
    forecast_horizon = forecast_horizon,
    calibration_time = calibration_time
  )

  lab_site_site_spine_od <- get_lab_site_site_spine(
    input_ww_data = recent_input_ww_data_for_stan
  )

  site_subpop_spine_od <- get_site_subpop_spine(
    input_ww_data = recent_input_ww_data_for_stan,
    input_count_data = input_count_data
  )

  lab_site_subpop_spine_od <- get_lab_site_subpop_spine(
    lab_site_site_spine = lab_site_site_spine,
    site_subpop_spine = site_subpop_spine_od
  )


  expect_error(get_stan_data(
    input_count_data,
    recent_input_ww_data_for_stan,
    date_time_spine,
    lab_site_site_spine_od,
    site_subpop_spine_od,
    lab_site_subpop_spine_od,
    last_count_data_date,
    first_count_data_date,
    forecast_date,
    forecast_horizon,
    calibration_time,
    generation_interval,
    inf_to_count_delay,
    infection_feedback_pmf,
    params,
    include_ww
  ))

  # Make wastewater data outside of scope of admissions data
  withr::with_seed(123, {
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
  })

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
  date_time_spine <- get_date_time_spine(
    forecast_date = forecast_date,
    input_count_data = input_count_data,
    last_count_data_date = last_count_data_date,
    forecast_horizon = forecast_horizon,
    calibration_time = calibration_time
  )
  lab_site_site_spine_old <- get_lab_site_site_spine(
    input_ww_data = old_input_ww_data_for_stan
  )

  site_subpop_spine_old <- get_site_subpop_spine(
    input_ww_data = old_input_ww_data_for_stan,
    input_count_data = input_count_data
  )

  lab_site_subpop_spine_old <- get_lab_site_subpop_spine(
    lab_site_site_spine = lab_site_site_spine_old,
    site_subpop_spine = site_subpop_spine_old
  )


  expect_error(get_stan_data(
    input_count_data,
    old_input_ww_data,
    date_time_spine,
    lab_site_site_spine_od,
    site_subpop_spine_od,
    lab_site_subpop_spine_od,
    last_count_data_date,
    first_count_data_date,
    forecast_date,
    forecast_horizon,
    calibration_time,
    generation_interval,
    inf_to_count_delay,
    infection_feedback_pmf,
    params,
    include_ww
  ))
})

test_that("Test that pmf check works as expected", {
  expect_warning(get_stan_data(
    input_count_data,
    input_ww_data,
    date_time_spine,
    lab_site_site_spine,
    site_subpop_spine,
    lab_site_subpop_spine,
    last_count_data_date,
    first_count_data_date,
    forecast_date,
    forecast_horizon,
    calibration_time,
    generation_interval = to_simplex(rep(1, 100)),
    inf_to_count_delay,
    infection_feedback_pmf,
    params,
    include_ww
  ))

  expect_warning(get_stan_data(
    input_count_data,
    input_ww_data,
    date_time_spine,
    lab_site_site_spine,
    site_subpop_spine,
    lab_site_subpop_spine,
    last_count_data_date,
    first_count_data_date,
    forecast_date,
    forecast_horizon,
    calibration_time,
    generation_interval,
    inf_to_count_delay = to_simplex(rep(1, 100)),
    infection_feedback_pmf,
    params,
    include_ww
  ))

  expect_warning(get_stan_data(
    input_count_data,
    input_ww_data,
    date_time_spine,
    lab_site_site_spine,
    site_subpop_spine,
    lab_site_subpop_spine,
    last_count_data_date,
    first_count_data_date,
    forecast_date,
    forecast_horizon,
    calibration_time,
    generation_interval,
    inf_to_count_delay,
    infection_feedback_pmf = to_simplex(rep(1, 100)),
    params,
    include_ww
  ))

  expect_error(get_stan_data(
    input_count_data,
    input_ww_data,
    date_time_spine,
    lab_site_site_spine,
    site_subpop_spine,
    lab_site_subpop_spine,
    last_count_data_date,
    first_count_data_date,
    forecast_date,
    forecast_horizon,
    calibration_time,
    generation_interval,
    inf_to_count_delay,
    infection_feedback_pmf = c(0.5, 0.4, 0.2),
    params,
    include_ww
  ))
})
