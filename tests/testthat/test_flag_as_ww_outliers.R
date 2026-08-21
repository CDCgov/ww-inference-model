# Create a dummy dataset that mimics the structure expected by the function
dummy_data <- tibble::tibble(
  site = rep("Site1", 10),
  lab = rep("Lab1", 10),
  lab_site_index = rep(1, 10),
  date = as.Date("2021-01-01") + 0:9,
  log_genome_copies_per_ml = log(c(
    100,
    150,
    100,
    200,
    270,
    200,
    NA,
    400,
    20,
    600
  )),
  below_lod = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
)


test_that("function returns a dataframe with correct columns", {
  result <- flag_ww_outliers(dummy_data)

  testthat::expect_true("flag_as_ww_outlier" %in% names(result))
  testthat::expect_true("exclude" %in% names(result))
})

test_that("function flags outliers correctly", {
  # Modify dummy_data to create an outlier scenario
  dummy_data <- tibble::tibble(
    site = rep("Site1", 12),
    lab = rep("Lab1", 12),
    lab_site_index = rep(1, 12),
    date = as.Date("2021-01-01") + 0:11,
    log_genome_copies_per_ml = log(c(
      100,
      120,
      100,
      110,
      115,
      130,
      110,
      200,
      NA,
      100,
      20,
      500000
    )),
    below_lod = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
  )

  result <- flag_ww_outliers(dummy_data, log_conc_threshold = 1)

  # Check if the known outlier is flagged correctly
  testthat::expect_true(sum(result$flag_as_ww_outlier) > 0)
  # Check that this hasn't yet been labeled for exclusion
  testthat::expect_true(sum(result$exclude) == 0)
})

test_that("function does not flag non-outliers", {
  # Modify dummy_data to have no outliers
  dummy_data <- tibble::tibble(
    site = rep("Site1", 12),
    lab = rep("Lab1", 12),
    lab_site_index = rep(1, 12),
    date = as.Date("2021-01-01") + 0:11,
    log_genome_copies_per_ml = log(c(
      100,
      120,
      100,
      110,
      115,
      130,
      110,
      200,
      NA,
      100,
      20,
      150
    )),
    below_lod = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
  )

  result <- flag_ww_outliers(dummy_data)

  # Check that there is no outlier to be flagged
  testthat::expect_true(sum(result$flag_as_ww_outlier) == 0)
  # Check that this hasn't yet been labeled for exclusion
  testthat::expect_true(sum(result$exclude) == 0)
})

test_that("function handles NA values appropriately", {
  # Include NAs in concentration column and check behavior
  dummy_data <- tibble::tibble(
    site = rep("Site1", 12),
    lab = rep("Lab1", 12),
    lab_site_index = rep(1, 12),
    date = as.Date("2021-01-01") + 0:11,
    log_genome_copies_per_ml = c(
      NA,
      120,
      100,
      110,
      NA,
      130,
      110,
      200,
      NA,
      100,
      20,
      150
    ),
    below_lod = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
  )

  result <- flag_ww_outliers(dummy_data)

  # Ensure NAs are handled according to specifications (not flagged as outliers)
  testthat::expect_true(sum(result$flag_as_ww_outlier) == 0)
  # Check that this hasn't yet been labeled for exclusion
  testthat::expect_true(sum(result$exclude) == 0)
})

test_that("rho_threshold and log_conc threshold parameters works as expected", {
  dummy_data <- tibble::tibble(
    site = rep("Site1", 12),
    lab = rep("Lab1", 12),
    lab_site_index = rep(1, 12),
    date = as.Date("2021-01-01") + 0:11,
    log_genome_copies_per_ml = c(
      100,
      120,
      100,
      110,
      115,
      1000,
      110,
      100,
      NA,
      100,
      20,
      100
    ),
    below_lod = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
  )

  # With a high rho threshold, we definitely won't flag the last result
  high_rho_threshold_result <- flag_ww_outliers(
    dummy_data,
    rho_threshold = 5,
    log_conc_threshold = 5
  )
  # With a low rho threshold we should flag the last point
  low_rho_threshold_result <- flag_ww_outliers(
    dummy_data,
    rho_threshold = 0.5,
    log_conc_threshold = 0.5
  )

  # Expect fewer outliers with higher threshold and more with lower threshold
  testthat::expect_true(
    sum(
      low_rho_threshold_result$flag_as_ww_outlier
    ) >
      sum(
        high_rho_threshold_result$flag_as_ww_outlier
      )
  )
})


test_that("exclude column is set to zero by default", {
  result <- flag_ww_outliers(dummy_data)

  expect_true(all(result$exclude == 0))
})

# The rho rule is meant to catch a spike: a steep move followed immediately by a
# steep move in the opposite direction. A steep move followed by another in the
# same direction is a real trend and should be left alone.
make_ww_series <- function(log_conc) {
  n <- length(log_conc)
  tibble::tibble(
    site = rep("Site1", n),
    lab = rep("Lab1", n),
    lab_site_index = rep(1, n),
    date = as.Date("2021-01-01") + seq_len(n) - 1,
    log_genome_copies_per_ml = log_conc,
    below_lod = rep(0, n)
  )
}

test_that("a sustained trend is not flagged as a rho outlier", {
  rising <- make_ww_series(c(1, 1, 1, 1, 6, 11, 11, 11))
  falling <- make_ww_series(c(11, 11, 11, 11, 6, 1, 1, 1))

  for (rho_threshold in c(0.8, 1, 1.2)) {
    rising_result <- flag_ww_outliers(
      rising,
      rho_threshold = rho_threshold,
      log_conc_threshold = 100
    )
    falling_result <- flag_ww_outliers(
      falling,
      rho_threshold = rho_threshold,
      log_conc_threshold = 100
    )

    testthat::expect_equal(sum(rising_result$flag_as_ww_outlier), 0)
    testthat::expect_equal(sum(falling_result$flag_as_ww_outlier), 0)
  }
})

test_that("a spike is still flagged as a rho outlier", {
  spike <- make_ww_series(c(1, 1, 1, 1, 9, 1, 1, 1))

  for (rho_threshold in c(0.8, 1, 1.2)) {
    result <- flag_ww_outliers(
      spike,
      rho_threshold = rho_threshold,
      log_conc_threshold = 100
    )

    testthat::expect_equal(sum(result$flag_as_ww_outlier), 1)
  }
})
