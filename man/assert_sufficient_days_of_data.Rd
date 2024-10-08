% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkers.R
\name{assert_sufficient_days_of_data}
\alias{assert_sufficient_days_of_data}
\title{Assert that the vector of dates spans at least the specified
calibration time}
\usage{
assert_sufficient_days_of_data(
  date_vector,
  data_name,
  calibration_time,
  call = rlang::caller_env(),
  add_err_msg = ""
)
}
\arguments{
\item{date_vector}{the vector of dates to check, must be of Date type}

\item{data_name}{What data correspond to the dates in \code{date_vector}.
Used to make the error message informative (e.g.
"hospital admissions data")}

\item{calibration_time}{integer indicating the number of days that
the dates must span}

\item{call}{Calling environment to be passed to \code{\link[cli:cli_abort]{cli::cli_abort()}} for
traceback.}

\item{add_err_msg}{add_err_msg string containing an additional error message,
default is the empty string (\code{""})}
}
\value{
NULL invisible
}
\description{
Assert that the vector of dates spans at least the specified
calibration time
}
