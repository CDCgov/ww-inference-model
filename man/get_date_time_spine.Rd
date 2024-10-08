% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_stan_data.R
\name{get_date_time_spine}
\alias{get_date_time_spine}
\title{Get date time spine to map to model output}
\usage{
get_date_time_spine(
  forecast_date,
  input_count_data,
  last_count_data_date,
  calibration_time,
  forecast_horizon
)
}
\arguments{
\item{forecast_date}{a character string in ISO8601 format (YYYY-MM-DD)
indicating the date that the forecast is to be made.}

\item{input_count_data}{a dataframe of the count data to be passed
directly to stan, , must have the following columns: date, count, total_pop}

\item{last_count_data_date}{string indicating the date of the last observed
count data point in 1SO8601 format (YYYY-MM-DD)}

\item{calibration_time}{integer indicating the number of days to calibrate
the model for, default is \code{90}}

\item{forecast_horizon}{integer indicating the number of days, including the
forecast date, to produce forecasts for, default is \code{28}}
}
\value{
a tibble containing an integer for time mapped to the corresponding
date, for the entire calibration and forecast period
}
\description{
Get date time spine to map to model output
}
