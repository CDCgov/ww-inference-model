test_that("Test that the model runs on simulated data without wastewater.", {
  #######
  # run model briefly on the simulated data
  #######
  model_test_data_no_ww <- model_test_data
  model_test_data_no_ww$model_spec$include_ww <- 0
  model_test_data_no_ww$ww_data <- c()

  withr::with_seed(5, {
    fit <- do.call(
      wwinference::wwinference,
      model_test_data_no_ww
    )
  })
})
