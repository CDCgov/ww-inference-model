# Test data setup
ww_data <- tibble::tibble(
  date = lubridate::ymd(rep(c("2023-11-01", "2023-11-02"), 2)),
  site = c(rep(1, 2), rep(2, 2)),
  lab = c(1, 1, 1, 1),
  conc = c(345.2, 784.1, 401.5, 681.8),
  lod = c(20, 20, 15, 15),
  site_pop = c(rep(1e6, 2), rep(3e5, 2))
)

# Test that function returns a dataframe with correct columns
test_that("Function returns dataframe with correct columns", {
  processed <- preprocess_ww_data(ww_data,
    conc_col_name = "conc",
    lod_col_name = "lod"
  )

  expected_cols <- c(
    "date", "site", "lab", "genome_copies_per_ml",
    "lod", "lab_site_index", "site_index",
    "flag_as_ww_outlier", "lab_site_name",
    "forecast_date"
  )

  expect_true(all(expected_cols %in% names(processed)))
})
