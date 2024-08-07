test_that(
  "Test matrix normalization stan function and its purpose.",
  {
    model <- compiled_site_inf_model
    space_model_fxns <- spatial_fxns$functions


    sigma2_generalized <- 0.2
    dist_matrix <- as.matrix(
      dist(
        data.frame(
          x = c(850, 370, 480, 70),
          y = c(120, 750, 810, 960),
          diag = TRUE,
          upper = TRUE
        )
      )
    )
    n <- ncol(
      dist_matrix
    )
    corr_matrix <- space_model_fxns$exponential_decay_corr_func(
      dist_matrix = dist_matrix,
      phi = 25,
      l = 1
    )
    normalized_corr_matrix <- space_model_fxns$matrix_normalization(
      corr_matrix
    )

    c <- sigma2_generalized^(1 / n)
    sigma_matrix <- c * normalized_corr_matrix

    gen_var <- det(
      sigma_matrix
    )

    testthat::expect_equal(
      gen_var,
      sigma2_generalized
    )
  }
)
