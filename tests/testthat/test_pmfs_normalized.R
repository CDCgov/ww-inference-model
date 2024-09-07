test_that("Test that bundled PMFs in the package data sum to 1", {
  model <- compiled_site_inf_model

  shedding_pdf <- model$functions$get_vl_trajectory(
    tpeak = 5,
    viral_peak = 5,
    duration_shedding = 17,
    n = 100
  )

  testthat::expect_equal(sum(shedding_pdf), 1.0)


  default_spec <- wwinference::get_model_spec()

  generation_interval <- default_spec$generation_interval
  testthat::expect_equal(sum(generation_interval), 1.0)


  inf_to_count_delay <- default_spec$inf_to_count_delay
  testthat::expect_equal(sum(inf_to_count_delay), 1.0)

  inf_feedback <- default_spec$infection_feedback_pmf
  testthat::expect_equal(sum(inf_feedback), 1.0)
})
