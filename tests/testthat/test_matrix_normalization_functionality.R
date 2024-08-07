test_that(
  "Test matrix normalization stan function purpose.",
  {
    model <- compiled_site_inf_model
    space_model_fxns <- spatial_fxns$functions

    sigma2_generalized <- 0.2
    dist_matrix1 <- as.matrix(
      dist(
        data.frame(
          x = c(85, 37, 48, 7),
          y = c(12, 75, 81, 96),
          diag = TRUE,
          upper = TRUE
        )
      )
    )
    dist_matrix2 <- as.matrix(
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
      dist_matrix1
    )
    c <- sigma2_generalized^(1 / n)

    corr_matrix1 <- space_model_fxns$exponential_decay_corr_func(
      dist_matrix = dist_matrix1,
      phi = 25,
      l = 1
    )
    normalized_corr_matrix1 <- space_model_fxns$matrix_normalization(
      corr_matrix1
    )
    sigma_matrix1 <- c * normalized_corr_matrix1

    corr_matrix2 <- space_model_fxns$exponential_decay_corr_func(
      dist_matrix = dist_matrix2,
      phi = 25,
      l = 1
    )
    normalized_corr_matrix2 <- space_model_fxns$matrix_normalization(
      corr_matrix2
    )
    sigma_matrix2 <- c * normalized_corr_matrix2


    testthat::expect_equal(
      det(sigma_matrix1),
      det(sigma_matrix2)
    )
  }
)
