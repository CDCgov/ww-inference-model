test_that(
  "Test Exponential Decay Corr Func in Stan agrees with R",
  {
    model <- compiled_site_inf_model
    space_model_fxns <- spatial_fxns

    corr_func_params <- list(
      dist_matrix = as.matrix(dist(
        data.frame(
          x = runif(5, 0, 100),
          y = runif(5, 0, 100)
        ),
        diag = TRUE,
        upper = TRUE
      )),
      phi = 25,
      l = 1
    )
    exp_decay_corr_func_stan <-
      space_model_fxns$functions$exponential_decay_corr_func(
        corr_func_params$dist_matrix,
        corr_func_params$phi,
        corr_func_params$l
      )
    exp_decay_corr_func_r <- exponential_decay_corr_func_r(
      corr_func_params
    )

    testthat::expect_equal(
      exp_decay_corr_func_stan,
      exp_decay_corr_func_r,
      ignore_attr = TRUE
    )
  }
)
