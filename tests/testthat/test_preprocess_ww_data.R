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
    "date", "site", "lab", "log_genome_copies_per_ml",
    "log_lod", "site_pop", "lab_site_index", "site_index",
    "flag_as_ww_outlier", "lab_site_name", "below_lod"
  )

  checkmate::expect_names(names(processed), must.include = expected_cols)
})

# Test that concentration column is renamed correctly
test_that("Concentration column is renamed correctly", {
  processed <- preprocess_ww_data(ww_data,
    conc_col_name = "conc",
    lod_col_name = "log_lod"
  )
  checkmate::expect_names(
    names(processed),
    must.include = "log_genome_copies_per_ml",
    disjunct.from = "conc"
  )
})

# Test that LOD column is renamed correctly
test_that("LOD column is renamed correctly", {
  ww_data_test <- ww_data |>
    dplyr::rename(LOD = lod)
  processed <- preprocess_ww_data(ww_data_test,
    conc_col_name = "conc",
    lod_col_name = "LOD"
  )

  checkmate::expect_names(names(processed),
    must.include = "log_lod",
    disjunct.from = "LOD"
  )
})

# Test that lab_site_index and site_index are created correctly
test_that("lab_site_index and site_index are created correctly", {
  processed <- preprocess_ww_data(ww_data,
    conc_col_name = "conc",
    lod_col_name = "lod"
  )

  expect_true(all(!is.na(processed$lab_site_index)))
  expect_true(all(!is.na(processed$site_index)))

  # Check for unique indices for each lab-site combination
  expect_equal(length(unique(processed$lab_site_index)), 2)

  # Check for unique indices for each site
  expect_equal(length(unique(processed$site_index)), 2)

  # More complex example where each sample is processed in a different lab,
  # so we will need 4 unique site identifiers

  test_ww_data <- tibble::tibble(
    date = lubridate::ymd(rep(c("2023-11-01", "2023-11-02"), 2)),
    site = c(rep(1, 2), 2, 1),
    lab = c(1, 2, 3, 4),
    conc = c(345.2, 784.1, 401.5, 681.8),
    lod = c(20, 20, 15, 15),
    site_pop = c(rep(1e6, 2), rep(3e5, 2))
  )

  processed <- preprocess_ww_data(test_ww_data,
    conc_col_name = "conc",
    lod_col_name = "lod"
  )

  expect_true(all(!is.na(processed$lab_site_index)))
  expect_true(all(!is.na(processed$site_index)))

  # Check for unique indices for each lab-site combination
  expect_equal(length(unique(processed$lab_site_index)), 4)

  # Check for unique indices for each site
  expect_equal(length(unique(processed$site_index)), 2)


  # Check to make sure that lab and site indices get created correctly
  # even if labs and sites are characters
  test_ww_data <- tibble::tibble(
    date = lubridate::ymd(rep(c("2023-11-01", "2023-11-02"), 2)),
    site = c(rep("site1", 2), rep("site2", 2)),
    lab = rep("lab1", 4),
    conc = c(345.2, 784.1, 401.5, 681.8),
    lod = c(20, 20, 15, 15),
    site_pop = c(rep(1e6, 2), rep(3e5, 2))
  )

  processed <- preprocess_ww_data(test_ww_data,
    conc_col_name = "conc",
    lod_col_name = "lod"
  )

  expect_true(all(!is.na(processed$lab_site_index)))
  expect_true(all(!is.na(processed$site_index)))

  # Check for unique indices for each lab-site combination
  expect_equal(length(unique(processed$lab_site_index)), 2)

  # Check for unique indices for each site
  expect_equal(length(unique(processed$site_index)), 2)

  # Check that the correct integers are generated for indices
  expected_lab_site_indices <- c(1, 1, 2, 2)
  expected_site_indices <- c(1, 1, 2, 2)

  expect_equal(processed$lab_site_index, expected_lab_site_indices)
  expect_equal(processed$site_index, expected_site_indices)
})

# Test that below_lod flag is set correctly
test_that("below_lod flag is set correctly", {
  processed <- preprocess_ww_data(ww_data,
    conc_col_name = "conc",
    lod_col_name = "lod"
  )

  # Check if below_lod is a binary indicator (0 or 1)
  expect_true(all(processed$below_lod %in% c(0, 1)))

  # Verify correctness of below_lod values based on test data
  expected_below_lod <- c(0, 0, 0, 0) # none of the concentrations are
  # below LOD in test data
  expect_equal(processed$below_lod, expected_below_lod)

  # Make something below the LOD and ensure it works correctly
  ww_data_test <- tibble::tibble(
    date = lubridate::ymd(rep(c("2023-11-01", "2023-11-02"), 2)),
    site = c(rep(1, 2), rep(2, 2)),
    lab = c(1, 1, 1, 1),
    conc = c(10, 784.1, 401.5, 681.8),
    lod = c(20, 20, 15, 15),
    site_pop = c(rep(1e6, 2), rep(3e5, 2))
  )

  processed <- preprocess_ww_data(ww_data_test,
    conc_col_name = "conc",
    lod_col_name = "lod"
  )

  # Check if below_lod is a binary indicator (0 or 1)
  expect_true(all(processed$below_lod %in% c(0, 1)))

  # Verify correctness of below_lod values based on test data
  expected_below_lod <- c(1, 0, 0, 0) # none of the concentrations are
  # below LOD in test data
  expect_equal(processed$below_lod, expected_below_lod)
})

# Test that lab_site_index and site_index are created correctly
test_that("lab_site_index and site_index are created correctly", {
  processed <- preprocess_ww_data(ww_data,
    conc_col_name = "conc",
    lod_col_name = "lod"
  )

  expect_true(all(!is.na(processed$lab_site_index)))
  expect_true(all(!is.na(processed$site_index)))

  # Check for unique indices for each lab-site combination
  expect_equal(length(unique(processed$lab_site_index)), 2)

  # Check for unique indices for each site
  expect_equal(length(unique(processed$site_index)), 2)
})

# Test that lab_site_name is constructed properly
test_that("lab_site_name is constructed properly", {
  processed <- preprocess_ww_data(ww_data,
    conc_col_name = "conc",
    lod_col_name = "lod"
  )

  expected_lab_site_names <- c(
    "Site: 1, Lab: 1", "Site: 1, Lab: 1",
    "Site: 2, Lab: 1", "Site: 2, Lab: 1"
  )

  expect_equal(processed$lab_site_name, expected_lab_site_names)
})

# Test that the function handles empty dataframes correctly
test_that("Function handles empty dataframes with an error", {
  empty_ww_data <- ww_data[FALSE, ]

  expect_error(preprocess_ww_data(empty_ww_data,
    conc_col_name = "conc",
    lod_col_name = "lod"
  ))
})

# Test that the function flags outliers correctly
# (assuming flag_ww_outliers works as expected)
test_that("Function flags outliers correctly", {
  # Add a potential outlier value for testing
  outlier_ww_data <- tibble::tibble(
    date = seq(
      from = lubridate::ymd("2023-11-01"),
      to = lubridate::ymd("2023-11-20"),
      by = "days"
    ),
    site = rep(1, 20),
    lab = rep(1, 20),
    conc = c(
      345.2, 335.3, 345.2, 335.3, 345.2, 335.3,
      345.2, 335.3, 334,
      345.2, 335.3, 345.2, 335.3, 345.2, 335.3,
      345.2, 335.3, 340, 334, 100000
    ),
    lod = rep(20, 20),
    site_pop = rep(1e6, 20)
  )

  processed <- preprocess_ww_data(outlier_ww_data,
    conc_col_name = "conc",
    lod_col_name = "lod"
  )

  # Assuming flag_as_ww_outlier is binary (0 or 1) and
  # the last entry is an outlier
  expect_true(all(processed$flag_as_ww_outlier %in% c(0, 1)))

  # Check if the last entry is flagged as an outlier
  expect_equal(tail(processed$flag_as_ww_outlier, n = 1), c(1))
})

# Test that all rows are preserved after preprocessing
test_that("All rows are preserved after preprocessing", {
  processed <- preprocess_ww_data(ww_data,
    conc_col_name = "conc",
    lod_col_name = "lod"
  )

  # The number of rows should remain unchanged after processing
  expect_equal(nrow(processed), nrow(ww_data))
})

# Test that no NA values are introduced during preprocessing
# (assuming input has no NAs)
test_that("No NA values are introduced during preprocessing", {
  processed <- preprocess_ww_data(ww_data,
    conc_col_name = "conc",
    lod_col_name = "lod"
  )

  # Check for any new NA values introduced in any column
  expect_true(all(!is.na(processed)))
})

# Test that the function can handle LOD values equal to concentration values
# (edge case)
test_that("Function handles LOD values equal to concentration values", {
  edge_case_ww_data <- ww_data |>
    dplyr::mutate(conc = lod) # Set concentration equal to LOD,
  # we expect this should get flagged as below LOD

  processed_edge_case <- preprocess_ww_data(edge_case_ww_data,
    conc_col_name = "conc",
    lod_col_name = "lod"
  )

  # Check if below_lod is set to 1 when concentration equals LOD
  expect_equal(processed_edge_case$below_lod, rep(1, nrow(edge_case_ww_data)))
})
