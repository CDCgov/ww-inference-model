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
