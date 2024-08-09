# Test data setup
data <- tibble::tibble(
  date = lubridate::ymd(c("2023-10-01", "2023-10-02")),
  genome_copies_per_mL = c(300, 3e6),
  flag_as_ww_outlier = c(0, 1),
  exclude = c(0, 0)
)

# Test that function returns a dataframe with an additional 'exclude' column
test_that("Function returns dataframe with 'exclude' column", {
  processed <- indicate_ww_exclusions(data,
    outlier_col_name = "flag_as_ww_outlier",
    remove_outliers = TRUE
  )

  expect_true("exclude" %in% names(processed))
})

# Test that outliers are correctly marked for exclusion when
# remove_outliers is TRUE
test_that("Outliers are marked for exclusion when remove_outliers is TRUE", {
  processed <- indicate_ww_exclusions(data,
    outlier_col_name = "flag_as_ww_outlier",
    remove_outliers = TRUE
  )

  expected_exclude <- c(0, 1) # Second row should be marked as an outlier

  expect_equal(processed$exclude, expected_exclude)
})

# Test that no rows are marked for exclusion when remove_outliers is FALSE
test_that("No rows are marked for exclusion when remove_outliers is FALSE", {
  processed <- indicate_ww_exclusions(data,
    outlier_col_name = "flag_as_ww_outlier",
    remove_outliers = FALSE
  )

  expected_exclude <- c(0, 0) # No rows should be excluded

  expect_equal(processed$exclude, expected_exclude)
})
