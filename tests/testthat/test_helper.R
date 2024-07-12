test_that("Make sure we can find and load files we need for other tests.", {
  testthat::expect_true(
    exists("toy_stan_data")
  )

  testthat::expect_true(
    exists("toy_stan_fit_last_draw")
  )


  # Compiled model object should exist in the workspace, with functions exposed
  testthat::expect_true(
    exists("compiled_site_inf_model")
  )

  testthat::expect_true(
    "CmdStanModel" %in% class(compiled_site_inf_model)
  )

  testthat::expect_no_error(
    compiled_site_inf_model$functions$convert_to_logmean(1.0, 1.0)
  )
})
