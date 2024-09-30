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
