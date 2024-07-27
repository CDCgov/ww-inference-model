#' Check that all dates in dataframe passed in are before a specified date
#'
#' @param df dataframe with `date` column
#' @param max_date string indicating the maximum date in ISO8601 convention
#' e.g. YYYY-MM-DD
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#'
#' @return NULL, invisibly
check_date <- function(df, max_date, call = rlang::caller_env()) {
  if (max(df$date) > max_date) {
    cli::cli_abort(
      c(
        "The data passed in has observations beyond the specified",
        "maximum date. Either this is the incorrect vintaged",
        "data, or the data needs to be filtered to only contain",
        "observations before the maximum date"
      ),
      call = call,
      class = "wwinference_input_data_error"
    )
  }
  invisible()
}


#' Check that R(t) specified for generating simulated data is of sufficient
#' length
#'
#' @param r_in_weeks a vector indicating the R(t) in weeks
#' @param ot integer indicating the observed time: length of hospital admissions
#'  calibration time in days
#' @param nt integer indicating the nowcast time: length of time between last
#' hospital admissions date and forecast date in days
#' @param forecast_horizon integer indicating the duration of the forecast in
#' days e.g. 28 days
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#'
#' @return NULL, invisible
check_rt_length <- function(r_in_weeks,
                            ot,
                            nt,
                            forecast_horizon,
                            call = rlang::caller_env()) {
  if (length(r_in_weeks) < (ot + nt + forecast_horizon) / 7) {
    cli::cli_abort(
      c(
        "The weekly R(t) specifed isn't long enough to produce",
        "infections for the specified observed/calibration time (`ot`)",
        "nowcast period (`nt`), or forecast horizon (`forecast_horizon`)"
      ),
      call = call,
      class = "wwinference_fwd_sim_specification_error"
    )
  }
  invisible()
}

#' Check that the sum of the wastewater site populations don't exceed
#' the total population
#'
#' @param pop_size integer indicating the population size in the hypothetical
#' state
#' @param ww_pop_sites vector indicating the population size in the
#' catchment area in each of those sites
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#'  traceback.
#'
#' @return NULL, invisibly
check_ww_site_pops <- function(pop_size,
                               ww_pop_sites,
                               call = rlang::caller_env()) {
  if (sum(ww_pop_sites) > pop_size) {
    cli::cli_abort(
      c(
        "The sum of the specified wastewater site populations is greater than",
        "the total population. Check to make sure population sizes",
        "are specified correctly",
        call = call,
        class = "wwinference_input_data_error"
      )
    )
  }
  invisible()
}

#' Check that the specified site and lab indices line up
#'
#' @param site vector of integers indicating which site (WWTP) each separate
#' lab-site observation comes frm
#' @param lab ector of integers indicating which lab the lab-site observations
#' come from
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#'
#' @return NULL, invisibly
check_site_and_lab_indices <- function(site,
                                       lab,
                                       call = rlang::caller_env()) {
  if (length(site) != length(lab)) {
    cli::cli_abort(
      c(
        "The specified site and lab indices don't align. The two",
        "indices should uniquely define the site-lab combinations"
      ),
      call = call,
      class = "wwinference_fwd_sim_specification_error"
    )
  }
  invisible()
}


#' Check that all elements of a vector are non-negative
#'
#' @param x vector of arguments to check for negativity
#' @param arg string to print the name of the element your checking
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#'
#' @return NULL, invisibly
check_elements_non_neg <- function(x, arg = "x", call = rlang::caller_env()) {
  # Greater than or equal to 0 or is NA
  is_non_neg <- (x >= 0) | is.na(x)
  if (!all(is_non_neg)) {
    cli::cli_abort(
      c("{.arg {arg}} has negative elements",
        "!" = "All elements must be 0 or greater",
        "i" = "Elements {.val {which(!is_non_neg)}} are negative"
      ),
      class = "wwinference_input_data_error",
      call = call
    )
  }
  invisible()
}
