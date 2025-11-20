t_length <- 127
forecast_date <- "2024-01-01"
data <- tibble::tibble(
  date = seq(
    from = lubridate::ymd("2023-10-01"),
    to = lubridate::ymd("2023-10-01") + lubridate::days(t_length - 1),
    by = "days"
  ),
  observed_value = sample(10:25, t_length, replace = TRUE)
)

draws <- tibble::tibble()
for (i in 1:100) {
  draws_i <- data |>
    dplyr::mutate(
      pred_value = observed_value +
        runif(t_length, min = -10, max = 10),
      draw = i
    )
  draws <- dplyr::bind_rows(draws, draws_i)
}

test_draws <- draws |>
  dplyr::mutate(
    observed_value = ifelse(date < forecast_date, observed_value, NA)
  )

test_eval_data <- data |>
  dplyr::rename("daily_hosp_admits_eval" = observed_value)


test_that("Test there is no error with eval data", {
  expect_no_error(
    get_plot_forecasted_counts(
      draws = test_draws,
      forecast_date = forecast_date,
      count_data_eval = test_eval_data,
      count_data_eval_col_name = "daily_hosp_admits_eval"
    )
  )
})


test_that("Test there is no error without eval data", {
  expect_no_error(
    get_plot_forecasted_counts(
      draws = test_draws,
      forecast_date = forecast_date
    )
  )
})
