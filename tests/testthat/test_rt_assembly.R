test_that(paste0(
  "Test that assembling an unadjusted R(t) ",
  "vector using R (language) code gives the ",
  "same result as doing it via our custom Stan",
  "functions when those are loaded into R"
), {
  model <- compiled_site_inf_model

  weeks_to_days <- get_ind_m(
    168,
    24
  )
  ndays <- dim(weeks_to_days)[1]
  nweeks <- dim(weeks_to_days)[2]

  ## Make sure we cover a wide range
  sigma <- 5
  ac <- 0.25

  withr::with_seed(5325, {
    std_normal <- rnorm((nweeks - 1))
    init_val <- rnorm(1, 0.5, 0.25)
  })

  ## Build the vector ourselves
  unadj_log_r_weeks_r <- diff_ar1_from_z_scores(
    init_val, ac, sigma, std_normal,
    stationary = FALSE
  )
  unadj_log_r_days_r <- rep(unadj_log_r_weeks_r,
    each = 7
  ) # In R, expand weekly to daily

  ## convert to linear scale and trim to size
  unadj_r_days_r <- exp(unadj_log_r_days_r)[1:ndays]


  ## Compute in the same way it is done in stan
  unadj_log_r_weeks_stan <- model$functions$diff_ar1(
    init_val,
    ac,
    sigma,
    std_normal,
    is_stat = FALSE
  ) |>
    as.numeric()

  expect_equal(
    unadj_log_r_weeks_r,
    unadj_log_r_weeks_stan
  )

  unadj_r_days_stan <- exp(
    weeks_to_days %*% unadj_log_r_weeks_stan
  ) |>
    as.numeric()

  expect_equal(
    unadj_r_days_r,
    unadj_r_days_stan
  )
})
