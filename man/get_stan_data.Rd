% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_stan_data.R
\name{get_stan_data}
\alias{get_stan_data}
\title{Get stan data for ww + hosp model}
\usage{
get_stan_data(
  input_count_data,
  input_ww_data,
  date_time_spine,
  lab_site_site_spine,
  site_subpop_spine,
  lab_site_subpop_spine,
  last_count_data_date,
  first_count_data_date,
  forecast_date,
  forecast_horizon,
  calibration_time,
  generation_interval,
  inf_to_count_delay,
  infection_feedback_pmf,
  params,
  include_ww,
  compute_likelihood = 1
)
}
\arguments{
\item{input_count_data}{tibble with the input count data needed for stan}

\item{input_ww_data}{tibble with the input wastewater data and indices
needed for stan}

\item{date_time_spine}{tibble mapping dates to time in days}

\item{lab_site_site_spine}{tibble mapping lab-sites to sites}

\item{site_subpop_spine}{tibble mapping sites to subpopulations}

\item{lab_site_subpop_spine}{tibble mapping lab-sites to subpopulations}

\item{last_count_data_date}{string indicating the date of the last data
point in the count dataset in ISO8601 convention e.g. YYYY-MM-DD}

\item{first_count_data_date}{string indicating the date of the first data
point in the count dataset in ISO8601 convention e.g. YYYY-MM-DD}

\item{forecast_date}{string indicating the forecast date in ISO8601
convention e.g. YYYY-MM-DD}

\item{forecast_horizon}{integer indicating the number of days to make a
forecast for}

\item{calibration_time}{integer indicating the max duration in days that
the model is calibrated to the count data for}

\item{generation_interval}{a vector with a zero-truncated normalized pmf of
the generation interval}

\item{inf_to_count_delay}{a vector with a normalized pmf of the delay from
infection to counts}

\item{infection_feedback_pmf}{a vector with a normalized pmf dictating the
delay of infection feedback}

\item{params}{a list mapping parameter names to their values}

\item{include_ww}{integer either 1 or 0 indicating whether to fit the
wastewater data or not, default is 1}

\item{compute_likelihood}{indicator variable telling stan whether or not to
compute the likelihood, default = \code{1}}
}
\value{
\code{stan_args}: named variables to pass to stan
}
\description{
Get stan data for ww + hosp model
}
\examples{
ww_data <- tibble::tibble(
  date = rep(seq(
    from = lubridate::ymd("2023-08-01"),
    to = lubridate::ymd("2023-11-01"),
    by = "weeks"
  ), 2),
  site = c(rep(1, 14), rep(2, 14)),
  lab = c(rep(1, 28)),
  conc = abs(rnorm(28, mean = 500, sd = 50)),
  lod = c(rep(20, 14), rep(15, 14)),
  site_pop = c(rep(2e5, 14), rep(4e5, 14))
)

ww_data_preprocessed <- preprocess_ww_data(ww_data,
  conc_col_name = "conc",
  lod_col_name = "lod"
)
input_ww_data <- indicate_ww_exclusions(ww_data_preprocessed)

hosp_data <- tibble::tibble(
  date = seq(
    from = lubridate::ymd("2023-07-01"),
    to = lubridate::ymd("2023-10-30"),
    by = "days"
  ),
  daily_admits = sample(5:70, 122, replace = TRUE),
  state_pop = rep(1e6, 122)
)

input_count_data <- preprocess_count_data(
  hosp_data,
  "daily_admits",
  "state_pop"
)

generation_interval <- to_simplex(c(0.01, 0.2, 0.3, 0.2, 0.1, 0.1, 0.01))
inf_to_count_delay <- to_simplex(c(
  rep(0.01, 12), rep(0.2, 4),
  rep(0.01, 10)
))
infection_feedback_pmf <- generation_interval

params <- get_params(
  system.file("extdata", "example_params.toml",
    package = "wwinference"
  )
)
forecast_date <- "2023-11-06"
calibration_time <- 90
forecast_horizon <- 28
include_ww <- 1
input_count_data_for_stan <- get_input_count_data_for_stan(
  input_count_data,
  calibration_time
)
last_count_data_date <- max(input_count_data_for_stan$date, na.rm = TRUE)
first_count_data_date <- min(input_count_data_for_stan$date, na.rm = TRUE)
input_ww_data_for_stan <- get_input_ww_data_for_stan(
  input_ww_data,
  first_count_data_date,
  last_count_data_date,
  calibration_time
)
date_time_spine <- get_date_time_spine(
  forecast_date = forecast_date,
  input_count_data = input_count_data_for_stan,
  last_count_data_date = last_count_data_date,
  forecast_horizon = forecast_horizon,
  calibration_time = calibration_time
)
lab_site_site_spine <- get_lab_site_site_spine(
  input_ww_data = input_ww_data_for_stan
)
site_subpop_spine <- get_site_subpop_spine(
  input_ww_data = input_ww_data_for_stan,
  input_count_data = input_count_data_for_stan
)
lab_site_subpop_spine <- get_lab_site_subpop_spine(
  lab_site_site_spine = lab_site_site_spine,
  site_subpop_spine
)


stan_data_list <- get_stan_data(
  input_count_data_for_stan,
  input_ww_data_for_stan,
  date_time_spine,
  lab_site_site_spine,
  site_subpop_spine,
  lab_site_subpop_spine,
  last_count_data_date,
  first_count_data_date,
  forecast_date,
  forecast_horizon,
  calibration_time,
  generation_interval,
  inf_to_count_delay,
  infection_feedback_pmf,
  params,
  include_ww
)
}
