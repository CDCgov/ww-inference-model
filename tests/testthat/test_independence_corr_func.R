test_that(
  "Test Independence Corr Func in Stan agrees with R",
  {
    model <- compiled_site_inf_model
    space_model_fxns <- spatial_fxns

    # corr func params 1
    corr_func_params <- list(
      num_sites = 5
    )
    ind_corr_func_stan <-
      space_model_fxns$functions$independence_corr_func(
        corr_func_params$num_sites
      )
    ind_corr_func_r <- independence_corr_func_r(
      corr_func_params
    )
    # both functions return correlation matrix using an
    # independent correlation function

    testthat::expect_equal(
      ind_corr_func_stan,
      ind_corr_func_r,
      ignore_attr = TRUE
    )


    # corr func params 2
    corr_func_params <- list(
      num_sites = 1
    )
    ind_corr_func_stan <-
      space_model_fxns$functions$independence_corr_func(
        corr_func_params$num_sites
      )
    ind_corr_func_r <- independence_corr_func_r(
      corr_func_params
    )
    # both functions return correlation matrix using an
    # independent correlation function

    testthat::expect_equal(
      ind_corr_func_stan,
      ind_corr_func_r,
      ignore_attr = TRUE
    )


    # corr func params 3
    corr_func_params <- list(
      num_sites = 500
    )
    ind_corr_func_stan <-
      space_model_fxns$functions$independence_corr_func(
        corr_func_params$num_sites
      )
    ind_corr_func_r <- independence_corr_func_r(
      corr_func_params
    )
    # both functions return correlation matrix using an
    # independent correlation function

    testthat::expect_equal(
      ind_corr_func_stan,
      ind_corr_func_r,
      ignore_attr = TRUE
    )
  }
)
