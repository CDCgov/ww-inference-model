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
  arg <- "ww_conc"
  if (!all(assertr::not_na(ww_conc))) {
    cli::cli_abort("NAs in {conc_col_name} column")
  }
  check_elements_non_neg(ww_conc, arg, call)
  checkmate::assert_vector(ww_conc)

  ww_lod <- ww_data |> dplyr::pull({
    lod_col_name
  })
  arg <- "ww_lod"
  if (!all(assertr::not_na(ww_lod))) {
    cli::cli_abort("NAs in LOD column")
  }
  check_elements_non_neg(ww_lod, arg, call)
  checkmate::assert_vector(ww_lod)

  # Wastewater date column should be of date type!
  ww_obs_dates <- ww_data$date
  arg <- "ww_obs_dates"
  if (!all(assertr::not_na(ww_obs_dates))) {
    cli::cli_abort("NAs in wastewater observation dates")
  }

  # Sites  either need to be integers or characters, not be missing, and be
  # non-negative
  site_labels <- ww_data$site
  arg <- "site_labels"
  check_int_or_char(site_labels, arg, call)
  if (!all(assertr::not_na(site_labels))) {
    cli::cli_abort("NAs in site column")
  }
  check_elements_non_neg(site_labels, arg, call)

  # Labs either need to be integers or characters, not be missing, and be
  # non-negative
  lab_labels <- ww_data$lab
  arg <- "lab_labels"
  check_int_or_char(lab_labels, arg, call)
  if (!all(assertr::not_na(lab_labels))) {
    cli::cli_abort("NAs in lab column")
  }
  check_elements_non_neg(lab_labels, arg, call)


  # Site populations should be integers, not be missing, and be
  # non-negative
  site_pops <- ww_data$site_pop
  arg <- "site_pops"
  checkmate::assert_integerish(site_pops)
  if (!all(assertr::not_na(site_pops))) {
    cli::cli_abort("NAs in site population size column")
  }
  check_elements_non_neg(site_pops, arg, call)


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
  check_elements_non_neg(counts, arg, call)
  checkmate::assert_vector(counts)
  checkmate::assert_integerish(counts)


  # Currently, the framework only supports a single population size for
  # an individual model fit. Therefore, check that there are not multiple
  # "global" population sizes being passed in.
  pop <- hosp_data |> dplyr::pull({
    pop_size_col_name
  })
  arg <- "global_pop"
  checkmate::check_integerish(pop)
  if (!all(assertr::not_na(pop))) {
    cli::cli_abort("NAs in {pop_size_col_name} column")
  }
  check_elements_non_neg(pop, arg, call)
  check_global_pop(pop, arg, call)


  # Date column should be of date type, for count data, there should only
  # be on observation per day
  count_dates <- hosp_data$date
  arg <- "count_obs_dates"
  checkmate::assert_date(count_dates)
  if (!all(assertr::is_uniq(count_dates))) {
    cli::cli_abort("count date dates are non-unique")
  }




  invisible()
}
