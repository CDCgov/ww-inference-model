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
  check_no_missingness(ww_conc, arg, call)
  check_elements_non_neg(ww_conc, arg, call)
  check_vector(ww_conc, arg, call)

  ww_lod <- ww_data |> dplyr::pull({
    lod_col_name
  })
  arg <- "ww_lod"
  check_no_missingness(ww_lod, arg, call)
  check_elements_non_neg(ww_lod, arg, call)
  check_vector(ww_lod, arg, call)

  # Wastewater date column should be of date type!
  ww_obs_dates <- ww_data$date
  arg <- "ww_obs_dates"
  check_date(ww_obs_dates, arg, call)
  check_no_missingness(ww_obs_dates, arg, call)

  # Sites  either need to be integers or characters, not be missing, and be
  # non-negative
  site_labels <- ww_data$site
  arg <- "site_labels"
  check_int_or_char(site_labels, arg, call)
  check_no_missingness(site_labels, arg, call)
  check_elements_non_neg(site_labels, arg, call)

  # Labs either need to be integers or characters, not be missing, and be
  # non-negative
  lab_labels <- ww_data$lab
  arg <- "lab_labels"
  check_int_or_char(lab_labels, arg, call)
  check_no_missingness(lab_labels, arg, call)
  check_elements_non_neg(lab_labels, arg, call)


  # Site populations should be integers, not be missing, and be
  # non-negative
  site_pops <- ww_data$site_pop
  arg <- "site_pops"
  check_int(site_pops, arg, call)
  check_no_missingness(site_pops, arg, call)
  check_elements_non_neg(site_pops, arg, call)


  invisible()
}
