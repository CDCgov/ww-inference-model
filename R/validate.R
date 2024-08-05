#' Validate user-provided wastewater concentration data
#'
#' @param ww_data tibble containing the input wastewater data
#' @param conc_col_name string indicating the name of the column containing
#' the concentration measurements in the wastewater data
#' @param lod_col_name string indicating the name of the column containing
#' the limit of detection for each measurement in the wastewater data
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#'
#' @return NULL, invisibly
validate_ww_conc_data <- function(ww_data,
                                  conc_col_name,
                                  lod_col_name,
                                  call = rlang::caller_env()) {
  ww_conc <- ww_data |> dplyr::pull({
    conc_col_name
  })
  arg <- conc_col_name
  assert_non_missingness(ww_conc, arg, call)
  assert_elements_non_neg(ww_conc, arg, call,
    add_err_msg = paste0(
      "Note that the model expects natural ",
      "scale concentration values, ",
      "which must be non-negative"
    )
  )
  checkmate::assert_vector(ww_conc)

  ww_lod <- ww_data |> dplyr::pull({
    lod_col_name
  })
  arg <- "lod_col_name"
  assert_non_missingness(ww_lod, arg, )
  assert_elements_non_neg(ww_lod, arg, call,
    add_err_msg = paste0(
      "Note that the model expects natural ",
      "scale LOD values, which must be ",
      "non-negative"
    )
  )
  checkmate::assert_vector(ww_lod)

  # Wastewater date column should be of date type!
  ww_obs_dates <- ww_data$date
  arg <- "ww_obs_dates"
  assert_non_missingness(ww_obs_dates, arg, call)
  checkmate::assert_date(ww_obs_dates)
  # check for multiple observations per day within a site and lab
  ww_data |>
    dplyr::group_by(site, lab) |>
    assert_no_repeated_elements()

  # Sites  either need to be integers or characters, not be missing, and be
  # non-negative
  site_labels <- ww_data$site
  arg <- "site_labels"
  assert_int_or_char(site_labels, arg, call)
  assert_non_missingness(site_labels, arg, call)

  # Labs either need to be integers or characters, not be missing, and be
  # non-negative
  lab_labels <- ww_data$lab
  arg <- "lab_labels"
  assert_int_or_char(lab_labels, arg, call)
  assert_non_missingness(lab_labels, arg, call)


  # Site populations should be integers, not be missing, and be
  # non-negative
  site_pops <- ww_data$site_pop
  arg <- "site_pops"
  checkmate::assert_integerish(site_pops)
  assert_non_missingness(site_pops, arg, call)
  assert_elements_non_neg(site_pops, arg, call)


  invisible()
}

#' Validate user-provided count data
#'
#' @param hosp_data tibble containing the input count data
#' @param count_col_name string indicating the name of the column containing
#' the count data
#' @param pop_size_col_name string indicating the name of the column containing
#' the population size of the count catchment area
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#'
#' @return NULL, invisibly
validate_count_data <- function(hosp_data,
                                count_col_name,
                                pop_size_col_name,
                                call = rlang::caller_env()) {
  # Count data should be non negative and a vector of integers
  counts <- hosp_data |> dplyr::pull({
    count_col_name
  })
  arg <- "counts"
  checkmate::assert_vector(counts)
  checkmate::assert_integerish(counts)
  assert_elements_non_neg(counts, arg, call)


  # Currently, the framework only supports a single population size for
  # an individual model fit. Therefore, check that there are not multiple
  # "global" population sizes being passed in.
  pop <- hosp_data |> dplyr::pull({
    pop_size_col_name
  })
  arg <- "global_pop"
  checkmate::check_integerish(pop)
  assert_elements_non_neg(pop)
  assert_non_missingness(pop, arg, call)
  add_err_msg <- paste0(
    "Multiple/time-varying count catchment area populations ",
    "are not currently supported. Check that data is from a ",
    "single location, and if so, consider replacing with an ",
    "average population size over the inference period"
  )
  assert_single_value(pop, arg, call, add_err_msg)


  # Date column should be of date type, for count data, there should only
  # be one observation per day
  count_dates <- hosp_data$date
  arg <- "count_obs_dates"
  checkmate::assert_date(count_dates)
  add_err_msg <- paste0(
    "Check that data is from a single location, and if so, ",
    "ensure that there are not multiple count data streams"
  )
  assert_no_repeated_elements(count_dates, arg, call, add_err_msg)





  invisible()
}
