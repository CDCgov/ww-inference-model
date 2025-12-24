#' Make daily data into weekly data
#' This is an internal function used to generate simulated weekly data
#'
#' @param hosp_data A tibble containing daily hospital admissions data,
#' expects the following columns: "date", and whatever is specified in
#' `count_col_name`
#' @param count_col_name A character string indicating the name of the column
#' with the daily counts, default is `"daily_hosp_admits"`
#' @param day_of_week_to_sum A character string indicating what day of the
#' week to assign the past 7 days of hospital admissions to. Must be full
#' weekday name with first letter in uppercase. Default is Saturday
#'
#' @return A dataframe with weekly hospital admissions data, assigned to
#' the day of the week to sum
make_weekly_data <- function(hosp_data,
                             count_col_name = "daily_hosp_admits",
                             day_of_week_to_sum = "Saturday") {
  hosp_data_w_wday <- hosp_data |>
    dplyr::mutate(
      day_of_week = lubridate::wday(.data$date,
        week_start = 1,
        label = TRUE,
        abbr = FALSE
      ),
      day_of_week_numeric = lubridate::wday(.data$date),
      weekly_hosp_admits = zoo::rollsum(.data[[count_col_name]],
        k = 7,
        na.pad = TRUE,
        align = "right"
      )
    ) |>
    dplyr::filter(day_of_week == {{ day_of_week_to_sum }})

  weekly_hosp_data <- hosp_data_w_wday |>
    dplyr::select(date, weekly_hosp_admits, state_pop)

  return(weekly_hosp_data)
}
