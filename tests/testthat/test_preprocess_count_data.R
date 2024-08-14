# Test data setup
hosp_data <- tibble::tibble(
  date = lubridate::ymd(c("2023-11-01", "2023-11-02")),
  daily_admits = c(10, 20),
  state_pop = c(1e6, 1e6)
)

# Test that function returns a dataframe with correct columns
test_that("Function returns dataframe with correct columns", {
  processed <- preprocess_count_data(hosp_data,
    count_col_name = "daily_admits",
    pop_size_col_name = "state_pop"
  )

  expected_cols <- c("date", "count", "total_pop")

  expect_true(all(expected_cols %in% names(processed)))
})

# Test that count column is renamed correctly
test_that("Count column is renamed correctly", {
  processed <- preprocess_count_data(hosp_data,
    count_col_name = "daily_admits",
    pop_size_col_name = "state_pop"
  )

  expect_false("daily_admits" %in% names(processed))
  expect_true("count" %in% names(processed))
})

# Test that population size column is renamed correctly
test_that("Population size column is renamed correctly", {
  processed <- preprocess_count_data(hosp_data,
    count_col_name = "daily_admits",
    pop_size_col_name = "state_pop"
  )

  expect_false("state_pop" %in% names(processed))
  expect_true("total_pop" %in% names(processed))
})

# Test data setup
hosp_data <- tibble::tibble(
  date = lubridate::ymd(c("2023-11-01", "2023-11-02")),
  daily_admits = c(10, 20),
  state_pop = c(1e6, 1e6)
)

# Test that function returns a dataframe with correct columns
test_that("Function returns dataframe with correct columns", {
  processed <- preprocess_count_data(hosp_data,
    count_col_name = "daily_admits",
    pop_size_col_name = "state_pop"
  )

  expected_cols <- c("date", "count", "total_pop")

  expect_true(all(expected_cols %in% names(processed)))
})

# Test that count column is renamed correctly
test_that("Count column is renamed correctly", {
  processed <- preprocess_count_data(hosp_data,
    count_col_name = "daily_admits",
    pop_size_col_name = "state_pop"
  )

  expect_false("daily_admits" %in% names(processed))
  expect_true("count" %in% names(processed))
})

# Test that population size column is renamed correctly
test_that("Population size column is renamed correctly", {
  processed <- preprocess_count_data(hosp_data,
    count_col_name = "daily_admits",
    pop_size_col_name = "state_pop"
  )

  expect_false("state_pop" %in% names(processed))
  expect_true("total_pop" %in% names(processed))
})

test_that("Function handles missing columns with an error", {
  incomplete_hosp_data <- hosp_data |> dplyr::select(-daily_admits)
  

  expect_error(preprocess_count_data(incomplete_hosp_data,
    count_col_name = "daily_admits",
    pop_size_col_name = "state_pop"
  ))
})

test_that("All rows are preserved after preprocessing", {
  processed <- preprocess_count_data(hosp_data,
    count_col_name = "daily_admits",
    pop_size_col_name = "state_pop"
  )

  # The number of rows should remain unchanged after processing
  expect_equal(nrow(processed), nrow(hosp_data))
})

# Test that no NA values are introduced during preprocessing
# (assuming input has no NAs)
test_that("No NA values are introduced during preprocessing", {
  processed <- preprocess_count_data(hosp_data,
    count_col_name = "daily_admits",
    pop_size_col_name = "state_pop"
  )

  # Check for any new NA values introduced in any column
  expect_true(all(complete.cases(processed)))
})

# Test that renaming columns works when using default parameter values
test_that("Renaming columns works with default parameters", {
  test_hosp_data <- tibble::tibble(
    date = lubridate::ymd(c("2023-11-01", "2023-11-02")),
    daily_hosp_admits = c(10, 20),
    state_pop = c(1e6, 1e6)
  )
  # Using default column names
  default_processed <- preprocess_count_data(test_hosp_data)

  expect_false(any(c(
    "daily_admits", "state_pop"
  ) %in% names(default_processed)))
  expect_true(all(c("count", "total_pop") %in% names(default_processed)))
})

# Test that the function handles empty dataframes correctly
test_that("Function handles empty dataframes with an error", {
  empty_hosp_data <- hosp_data[FALSE, ]

  expect_error(preprocess_count_data(empty_hosp_data,
    count_col_name = "daily_admits",
    pop_size_col_name = "state_pop"
  ))
})

# Test that the function can handle zero counts and population sizes
# without errors
test_that("Function handles zero counts and population sizes without errors", {
  zero_counts_hosp_data <- hosp_data |>
    dplyr::mutate(daily_admits = 0, state_pop = 0)

  processed_zero_counts <- preprocess_count_data(zero_counts_hosp_data,
    count_col_name = "daily_admits",
    pop_size_col_name = "state_pop"
  )

  # Check if zeros are preserved in count and total_pop columns after processing
  expect_equal(
    processed_zero_counts$count,
    rep(0, nrow(zero_counts_hosp_data))
  )
  expect_equal(
    processed_zero_counts$total_pop,
    rep(0, nrow(zero_counts_hosp_data))
  )
})

# Test that the function can handle NAs in count column
test_that("Function handles NAs in count column", {
  na_counts_hosp_data <- hosp_data |>
    dplyr::mutate(daily_admits = c(NA, daily_admits[2]))

  processed_na_counts <- preprocess_count_data(na_counts_hosp_data,
    count_col_name = "daily_admits",
    pop_size_col_name = "state_pop"
  )

  # Check if NAs are preserved in count columns after processing
  expect_true(any(is.na(processed_na_counts$count)))
})
