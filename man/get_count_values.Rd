% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_stan_data.R
\name{get_count_values}
\alias{get_count_values}
\title{Get count values}
\usage{
get_count_values(input_count_data, ot, ht, count_col_name = "count")
}
\arguments{
\item{input_count_data}{a dataframe with the input count data}

\item{ot}{integer indicating the total duration of time that the
model has available calibration data in days}

\item{ht}{integer indicating the number of days to produce count estimates
outside the calibration period (forecast + nowcast time) in days}

\item{count_col_name}{A string representing the name of the column in the
input_count_data that indicates the number of daily counts of the
epidemiological indicator, e.g. cases or hospital admissions,
default is \code{count}}
}
\value{
A list containing the necessary vectors of values that
the stan model requires:
counts: a vector of number of daily count observations
day_of_week: a vector indicating the day of the week of each of the dates
in the calibration and forecast period
}
\description{
Get count values
}
