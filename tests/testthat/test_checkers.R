test_that(
  "Test that check on maximum date errors as expected.",
  {
    date_vector <- seq(
      from = lubridate::ymd("2024-01-01"),
      to = lubridate::ymd("2024-01-05"),
      by = "days"
    )
    max_date <- lubridate::ymd("2024-01-05")

    expect_no_error(assert_no_dates_after_max(date_vector, max_date))

    max_date <- lubridate::ymd("2024-01-02")

    expect_error(assert_no_dates_after_max(date_vector, max_date))

    max_date <- "character"

    expect_error(assert_no_dates_after_max(date_vector, max_date))
  }
)

test_that(
  "Test that check on elements being non-negative works as expected.",
  {
    x <- c(1, 0, 7)
    expect_no_error(assert_elements_non_neg(x, "x"))

    x <- c(-9, 4, 5)
    expect_error(assert_elements_non_neg(x, "x"))
  }
)

test_that(
  "Test that check on non-missingness works as expected.",
  {
    x <- c(1, 0, 7)
    expect_no_error(assert_non_missingness(x, "x"))

    x <- c(-9, 4, NA)
    expect_error(assert_elements_non_neg(x, "x"))
  }
)

test_that(
  "Test that check on elements not repeating works as expected.",
  {
    x <- c(1, 0, 7)
    expect_no_error(assert_no_repeated_elements(x, "x"))

    x <- c(4, 4, 5)
    expect_error(assert_no_repeated_elements(x, "x"))

    x <- c("2024-03-03", "2024-03-03", "2024-03-04")
    expect_error(assert_no_repeated_elements(x, "x"))
  }
)

test_that(
  "Test that check on interger or characters works expected.",
  {
    x <- c(1, 0, 7)
    expect_no_error(assert_int_or_char(x, "x"))
    x <- c("site1", "site2")
    expect_no_error(assert_int_or_char(x, "x"))
    x <- c("site1", 4)
    expect_no_error(assert_int_or_char(x, "x"))


    x <- c(1.1, 3)
    expect_error(assert_int_or_char(x, "x"))
  }
)

test_that(
  "Test that check on wastewater data column names works as expected.",
  {
    x <- tibble::tibble(
      date = lubridate::ymd("2024-01-01"),
      site = 1,
      lab = 2,
      genome_copies_per_ml = 4,
      lod = 6,
      site_pop = 100,
      test_col = 4
    )
    conc_col_name <- "genome_copies_per_ml"
    lod_col_name <- "lod"
    expect_no_error(check_req_ww_cols_present(
      x,
      conc_col_name,
      lod_col_name
    ))

    # Pass wrong col name
    x <- tibble::tibble(
      date = lubridate::ymd("2024-01-01"),
      site = 1,
      lab = 2,
      genome_copies_per_ml = 4,
      lod = 6,
      site_pop = 100
    )
    conc_col_name <- "genome_copies_per_ml"
    lod_col_name <- "LOD"
    expect_error(assert_req_ww_columns_present(
      x,
      conc_col_name,
      lod_col_name
    ))

    # Wrong column name
    x <- tibble::tibble(
      date = lubridate::ymd("2024-01-01"),
      site = 1,
      lab = 2,
      genome_copies_per_ml = 4,
      lod = 6,
      site_pops = 100
    )
    conc_col_name <- "genome_copies_per_ml"
    lod_col_name <- "lod"
    expect_error(assert_req_ww_columns_present(
      x,
      conc_col_name,
      lod_col_name
    ))
  }
)

test_that(
  "Test that check on count data column names works as expected.",
  {
    x <- tibble::tibble(
      date = lubridate::ymd("2024-01-01"),
      hosp = 5,
      pop = 100
    )
    count_col_name <- "hosp"
    pop_size_col_name <- "pop"
    expect_no_error(check_req_count_cols_present(
      x,
      count_col_name,
      pop_size_col_name
    ))

    # Pass wrong column name
    x <- tibble::tibble(
      date = lubridate::ymd("2024-01-01"),
      hosp = 5,
      pop = 100
    )
    count_col_name <- "count"
    pop_size_col_name <- "pop"
    expect_error(check_req_hosp_columns_present(
      x,
      count_col_name,
      pop_size_col_name
    ))

    # Wrong column name
    x <- tibble::tibble(
      dates = lubridate::ymd("2024-01-01"),
      hosp = 5,
      pop = 100
    )
    count_col_name <- "hosp"
    pop_size_col_name <- "pop"
    expect_error(check_req_hosp_columns_present(
      x,
      count_col_name,
      pop_size_col_name
    ))
  }
)

test_that(
  "Test that check on only a single value being present works as expected.",
  {
    x <- c(1)
    expect_no_error(assert_single_value(x, "x"))

    x <- c(-9, 4, 5)
    expect_error(assert_single_value(x, "x"))
  }
)

test_that(
  "Test that check on non-empty dataframe works as expected.",
  {
    df_with_cols <- data.frame(Col1 = c(1, 2), Column2 = c(3, 4))
    expect_no_error(assert_df_not_empty(df_with_cols, "df_with_cols"))

    empty_df_with_cols <- data.frame(Col1 = numeric(), Column2 = character())
    expect_error(assert_df_not_empty(empty_df_with_cols, "empty_df_with_cols"))
  }
)

test_that(
  "Test that check on non-empty tibble works as expected.",
  {
    tibble_with_cols <- tibble::tibble(Col1 = c(1, 2), Column2 = c(3, 4))
    expect_no_error(assert_df_not_empty(tibble_with_cols, "tibble_with_cols"))

    empty_df_with_cols <- data.frame(Col1 = numeric(), Column2 = character())
    expect_error(assert_df_not_empty(empty_df_with_cols, "empty_df_with_cols"))
  }
)


test_that(
  "Test that check on daily dates works as expected.",
  {
    daily_dates <- c(
      lubridate::ymd("2023-01-01"),
      lubridate::ymd("2023-01-02")
    )
    expect_no_error(assert_daily_data(daily_dates))

    weekly_dates <- c(
      lubridate::ymd("2023-01-01"),
      lubridate::ymd("2023-01-08")
    )
    expect_error(assert_daily_data(weekly_dates))
  }
)

test_that(
  "Test that assert dates in range function works as expected.",
  {
    dates1 <- lubridate::ymd(c("2023-01-01", "2023-01-02"))
    dates2 <- lubridate::ymd(c("2023-01-01", "2023-01-04"))
    max_date <- "2023-01-05"
    expect_no_error(assert_dates_within_frame(
      dates1,
      dates2,
      max_date
    ))


    dates1 <- lubridate::ymd(c("2023-01-01", "2023-01-02"))
    dates2 <- lubridate::ymd(c("2023-01-03", "2023-01-04"))
    max_date <- "2023-01-05"
    expect_no_error(assert_dates_within_frame(
      dates1,
      dates2,
      max_date
    ))

    dates1 <- lubridate::ymd(c("2023-01-01", "2023-01-02"))
    dates2 <- lubridate::ymd(c("2024-01-03", "2024-01-04"))
    max_date <- "2023-01-05"
    expect_error(assert_dates_within_frame(
      dates1,
      dates2,
      max_date
    ))
  }
)
