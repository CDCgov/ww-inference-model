#' Check that all dates in dataframe passed in are before a specified date
#'
#' @param df dataframe with `date` column
#' @param max_date string indicating the maximum date in ISO8601 convention
#' e.g. YYYY-MM-DD
#' @param call Calling environment to be passed to the type checker
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


#' Check that all elements of a vector are non-negative
#'
#' @param x vector of arguments to check for negativity
#' @param arg string to print the name of the element your checking
#' @param call Calling environment to be passed to the type checker
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

#' Check that the input wastewater data contains all the required column names
#'
#' @param ww_data tibble containing the input wastewater data
#' @param conc_col_name string indicating the name of the column containing
#' the concentration measurements in the wastewater data
#' @param lod_col_name string indicating the name of the column containing
#' the concentration measurements in the wastewater data
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#'
#' @return NULL, invisibly
check_required_ww_inputs <- function(ww_data,
                                     conc_col_name,
                                     lod_col_name,
                                     call = rlang::caller_env()) {
  column_names <- colnames(ww_data)
  if (!"date" %in% column_names) {
    cli::cli_abort(
      c("`date` column missing from input wastewater data"),
      class = "wwinference_input_data_error",
      call = call
    )
  }

  if (!"site" %in% column_names) {
    cli::cli_abort(
      c("site column missing from input wastewater data"),
      class = "wwinference_input_data_error",
      call = call
    )
  }
  if (!"lab" %in% column_names) {
    cli::cli_abort(
      c("`lab` column missing from input wastewater data"),
      class = "wwinference_input_data_error",
      call = call
    )
  }
  if (!"site_pop" %in% column_names) {
    cli::cli_abort(
      c("`site_pop` column missing from input wastewater data"),
      class = "wwinference_input_data_error",
      call = call
    )
  }
  if (! !!conc_col_name %in% column_names) {
    cli::cli_abort(
      c("{conc_col_name} column missing from input wastewater data"),
      class = "wwinference_input_data_error",
      call = call
    )
  }
  if (! !!lod_col_name %in% column_names) {
    cli::cli_abort(
      c("{lod_col_name} column missing from input wastewater data"),
      class = "wwinference_input_data_error",
      call = call
    )
  }

  invisible()
}

#' Check there is no missignness in a particular vector
#'
#' @param x the vector to check
#' @param arg the name of the vector to check
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#'
#' @return NULL, invisibly
check_no_missingness <- function(x, arg = "x", call = rlang::caller_env()) {
  is_missing <- rlang::are_na(x)

  if (any(is_missing)) {
    cli::cli_abort(
      c("{.arg {arg}} has missing values",
        "i" = "Missing values are not supported in {.arg {arg}}",
        "!" = "Missing element(s) index: {.val {which(is_missing)}}"
      ),
      call = call,
      class = "wwinference_input_data_error"
    )
  }
}
