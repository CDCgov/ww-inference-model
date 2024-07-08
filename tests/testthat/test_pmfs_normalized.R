test_that("PMFs sum to 1", {
  model <- compiled_site_inf_model

  shedding_pdf <- model$functions$get_vl_trajectory(
    tpeak = 5,
    viral_peak = 5,
    duration_shedding = 17,
    n = 100
  )

  testthat::expect_equal(sum(shedding_pdf), 1.0)


  generation_interval <- toy_stan_data$generation_interval
  testthat::expect_equal(sum(generation_interval), 1.0)


  inf_to_count_delay <- toy_stan_data$inf_to_hosp
  testthat::expect_equal(sum(inf_to_count_delay), 1.0)

  inf_feedback <- toy_stan_data$infection_feedback_pmf
  testthat::expect_equal(sum(inf_feedback), 1.0)
})
