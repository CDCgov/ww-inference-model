test_that(
  "Test matrix normalization stan function.",
  {
    model <- compiled_site_inf_model
    space_model_fxns <- spatial_fxns$functions

    dist_matrix <- as.matrix(
      dist(
        data.frame(
          x = c(85, 37, 48, 7),
          y = c(12, 75, 81, 96),
          diag = TRUE,
          upper = TRUE
        )
      )
    )
    corr_matrix <- space_model_fxns$exponential_decay_corr_func(
      dist_matrix = dist_matrix,
      phi = 25,
      l = 1
    )
    normalized_corr_matrix <- space_model_fxns$matrix_normalization(
      corr_matrix
    )

    testthat::expect_equal(
      det(normalized_corr_matrix),
      1
    )
  }
)
