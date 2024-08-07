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
  }
)
