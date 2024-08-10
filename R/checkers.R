#' Check that all dates in dataframe passed in are before a specified date
#'
#' @description
#' This function is specifically meant to ensure that the data in the
#' `date_vector` specified does not contain dates after a given `max_date`. The
#' intended use-case for this is to ensure that one doesn't accidentally
#' pass in data that extends beyond the forecast date, as ideally the user
#' is providing vintaged "as of" datasets or at the very least is filtering the
#' data so that they are not including in their inference data that was
#' made available after the forecast was made.
#'
#'
#' @param date_vector vector of dates
#' @param max_date string indicating the maximum date in ISO8601 convention
#' e.g. YYYY-MM-DD
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#'
#' @return NULL, invisibly
assert_no_dates_after_max <- function(date_vector,
                                      max_date, call = rlang::caller_env()) {
  if (max(date_vector) > max_date) {
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


#' Assert that all elements of a vector are non-negative
#'
#' @param x vector of arguments to check for negativity
#' @param arg string to print the name of the element your checking
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#' @param add_err_msg string containing an additional error message,
#'  default is the empty string (`""`)
#'
#' @return NULL, invisibly
assert_elements_non_neg <- function(x, arg = "x",
                                    call = rlang::caller_env(),
                                    add_err_msg = "") {
  # Greater than or equal to 0 or is NA
  is_non_neg <- (x >= 0) | is.na(x)
  if (!all(is_non_neg)) {
    cli::cli_abort(
      c("{.arg {arg}} has negative elements.", add_err_msg,
        "!" = "All elements must be 0 or greater",
        "i" = "Elements {.val {which(!is_non_neg)}} are negative"
      ),
      class = "wwinference_input_data_error",
      call = call
    )
  }
  invisible()
}

#' Assert that there is no missignness in a particular vector
#'
#' @param x the vector to check
#' @param arg the name of the vector to check
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#' @param add_err_msg string containing an additional error message,
#' default is the empty string (`""`)
#'
#' @return NULL, invisibly
assert_non_missingness <- function(x, arg = "x",
                                   call = rlang::caller_env(),
                                   add_err_msg = "") {
  if (checkmate::anyMissing(x)) {
    cli::cli_abort(
      c("{.arg {arg}} has missing values", add_err_msg,
        "!" = "All elements of{.arg {arg}} should be present",
        "i" = "Element(s) {.val {which(is.na(x))}} are missing"
      ),
      class = "wwinference_input_data_error",
      call = call
    )
  }
  invisible()
}

#' Check that there are no repeated elements in the vector
#'
#' @description
#' This function  checks that the elements of a vector are
#' not repeated.
#'
#' @param x the vector to check
#' @param arg the name of the vector to check
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#' @param add_err_msg string containing an additional error message,
#' default is the empty string (`""`)
#'
#' @return NULL, invisibly
assert_no_repeated_elements <- function(x, arg = "x",
                                        call = rlang::caller_env(),
                                        add_err_msg = "") {
  duplicates <- duplicated(x)
  if (any(duplicates)) {
    cli::cli_abort(
      c("{.arg {arg}} has more than one element", add_err_msg,
        "i" = "Multiple {.arg {arg}} are not currently supported.",
        "!" = "Duplicate element(s) index: {.val {which(duplicates)}}"
      ),
      call = call,
      class = "wwinference_input_data_error"
    )
  }
  invisible()
}



#' Assert that a vector is either of a vector of integers or a vector of
#' characters
#' @description
#' This is for unique identifiers of groupings, which we will allow to
#' either be character or integers.
#'
#'
#' @param x Object with type checking
#' @param arg Name of the argument supplying the object
#' @param call Calling environment to be passed to [cli::cli_abort()]
#'
#' @return NULL, invisibly
assert_int_or_char <- function(x, arg = "x", call = rlang::caller_env()) {
  # Check if its a character, if it is, check passes. If not,
  # check if its an integer
  is_integerish <- checkmate::check_integerish(x)
  is_character <- checkmate::check_character(x)
  is_int_or_char <- isTRUE(is_integerish) || isTRUE(is_character)

  if (!is_int_or_char) {
    throw_type_error(
      object = x,
      arg_name = arg,
      expected_type = "integer or character",
      call = call
    )
  }
  invisible()
}


#' Throw an informative type error on a user-provided input
#'
#' Follows the guidance from [rlang::abort()] on applying a call and a class
#' in the error message. Used as a base in type-checkers to throw a properly
#' formatted error.
#'
#' @param object Object with incorrect type
#' @param arg_name Name of the argument corresponding to `object`
#' @param expected_type The type that the user should provide instead
#' @param call The calling environment to be reflected in the error message
#'
#' @return This function is called for its side effect of throwing an error. It
#'   should never return.
throw_type_error <- function(object,
                             arg_name,
                             expected_type,
                             call = rlang::caller_env()) {
  cli::cli_abort(
    c("{.arg {arg_name}} is {.obj_type_friendly {object}}",
      "i" = "Must be of type {.emph {expected_type}}"
    ),
    call = call,
    class = "wwinference_type_error"
  )
}


#' Check that the input wastewater data contains all the required column names
#'
#' @description
#' This function is intended to be used to check that the wastewater data that
#' gets passed into [preprocess_ww_data()] contains the required columns. If
#' it does not, we want to tell the user which column are missing. This will not
#' however, ensure that the elements of the column are of the right type,
#' or check that the values of them make sense.
#'
#'
#' @param ww_data tibble containing the input wastewater data
#' @param conc_col_name string indicating the name of the column containing
#' the concentration measurements in the wastewater data
#' @param lod_col_name string indicating the name of the column containing
#' the concentration measurements in the wastewater data
#' @param add_req_col_names vector of strings indicating the required wastewater
#' column names, the defaults are `c("date", "site", "lab", "site_pop")`
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#'
#' @return NULL, invisibly
check_req_ww_cols_present <- function(ww_data,
                                      conc_col_name,
                                      lod_col_name,
                                      add_req_col_names = c(
                                        "date", "site",
                                        "lab", "site_pop"
                                      ),
                                      call = rlang::caller_env()) {
  column_names <- colnames(ww_data)
  expected_col_names <- c(
    {
      conc_col_name
    },
    {
      lod_col_name
    },
    add_req_col_names
  )

  # This either returns TRUE or tells you whats missing.
  name_check_result <- checkmate::check_names(column_names,
    must.include = expected_col_names
  )
  if (!isTRUE(name_check_result)) {
    cli::cli_abort(
      message = c(
        "Required columns are missing from the wastewater data. ",
        autoescape_brackets(name_check_result)
      ),
      class = "wwinference_input_data_error",
      call = call
    )
  }

  invisible()
}

#' Check that the input count data contains all the required column names
#'
#' @description
#' This function is intended to be used to check that the count data that
#' gets passed into [preprocess_count_data()] contains the required columns. If
#' it does not, we want to tell the user which columns are missing. This will
#' not however, ensure that the elements of the column are of the right type,
#' or check that the values of them make sense.
#'
#'
#' @param count_data tibble containing the input count data
#' @param count_col_name string indicating the name of the column containing
#' the count data
#' @param pop_size_col_name string indicating the name of the column containing
#' the population size of the count catchment area
#' @param add_req_col_names vector of strings indicating the required count
#' data column names, the defaults is `"date"`
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#'
#' @return NULL, invisibly
check_req_count_cols_present <- function(count_data,
                                         count_col_name,
                                         pop_size_col_name,
                                         add_req_col_names = c("date"),
                                         call = rlang::caller_env()) {
  column_names <- colnames(count_data)
  expected_col_names <- c(
    {
      count_col_name
    },
    {
      pop_size_col_name
    },
    add_req_col_names
  )

  # This tells you whats missing
  check_colnames <- checkmate::check_names(column_names,
    must.include = expected_col_names
  )

  # This tells you from where it is missing
  if (!isTRUE(check_colnames)) {
    cli::cli_abort(
      c(
        "Required columns are missing from the input count data",
        autoescape_brackets(check_colnames)
      ),
      class = "wwinference_input_data_error",
      call = call
    )
  }

  invisible()
}

#' Assert that there is only a single value in a particular column
#'
#' @param x the vector to check
#' @param arg the name of the vector to check
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#' @param add_err_msg string containing an additional error message,
#' default is the empty string (`""`)
#'
#' @return NULL, invisibly
assert_single_value <- function(x, arg = "x",
                                call = rlang::caller_env(),
                                add_err_msg = "") {
  unique_elements <- unique(x)

  if (length(unique_elements) > 1) {
    cli::cli_abort(
      c(
        "{.arg {arg}} should have a single unique value;",
        "it has more than one. ",
        add_err_msg
      ),
      call = call,
      class = "wwinference_input_data_error"
    )
  }
  invisible()
}

#' Assert that the dataframe being passed to the function is not empty
#'
#' @param x the dataframe to check
#' @param arg the name of the dataframe to check
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#' @param add_err_msg add_err_msg string containing an additional error message,
#' default is the empty string (`""`)
#'
#' @return NULL invisible
assert_not_empty <- function(x,
                             arg = "x",
                             call = rlang::caller_env(),
                             add_err_msg = "") {
  nrows <- nrow(x)

  if (is.null(nrows)) {
    cli::cli_abort(
      c(
        "Expected something with rows,",
        "i.e. 2-d array or dataframe-like. ",
        add_err_msg
      ),
      call = call,
      class = "wwinference_input_data_error"
    )
  } else if (nrows < 1) {
    cli::cli_abort(c("{.arg {arg}} is empty", add_err_msg),
      call = call,
      class = "wwinference_input_data_error"
    )
  }
  invisible()
}


#' Assert that the vector of dates being passed in contains dates for each day
#'
#' @description
#' This function checks to make sure that the date vector being passed in
#' is complete for every day between the minimum and maximum dates. It can
#' have repeated values.
#'
#' @param dates the vector of dates to check, must be of Date type
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#' @param add_err_msg add_err_msg string containing an additional error message,
#' default is the empty string (`""`)
#'
#' @return NULL invisible
assert_daily_data <- function(dates,
                              call = rlang::caller_env(),
                              add_err_msg = "") {
  # Generate a sequence of dates from the minimum to the
  # maximum date in the dataset
  expected_dates <- seq.Date(
    from = min(dates),
    to = max(dates),
    by = "day"
  )


  if (!all(expected_dates %in% sorted_dates)) {
    cli::cli_abort(
      c(
        "Vector of dates does not contain dates for each day",
        add_err_msg
      ),
      call = call,
      class = "wwinference_input_data_error"
    )
  }
  invisible()
}

#' Assert that the vector of dates spans at least the specified
#' calibration time
#'
#' @param date_vector the vector of dates to check, must be of Date type
#' @param calibration_time integer indicating the number of days that
#' the dates must span
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#' @param add_err_msg add_err_msg string containing an additional error message,
#' default is the empty string (`""`)
#'
#' @return NULL invisible
assert_sufficient_days_of_data <- function(date_vector,
                                           calibration_time,
                                           call = rlang::caller_env(),
                                           add_err_msg = "") {
  # check that you have sufficient count data for the calibration time
  calibration_start <- max(date_vector,
    na.rm = TRUE
  ) - lubridate::days(calibration_time) + 1
  check_sufficient_data <- min(date_vector, na.rm = TRUE) <= calibration_start
  if (!check_sufficient_data) {
    cli::cli_abort(
      c(
        "Insufficient data for specified calibration time"
      ),
      call = call,
      class = "wwinference_specification_error"
    )
  }
  invisible()
}

#' Assert that two vectors of dates overlap
#'
#' @param dates1 the vector of dates to check, must of date type
#' @param dates2 the vector of dates to compare to, must be of date type
#' @param max_date the maximum date the testing dates should be
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#' @param add_err_msg add_err_msg string containing an additional error message,
#' default is the empty string (`""`)
#'
#' @return NULL invisible
assert_overlap_dates <- function(dates1,
                                 dates2,
                                 max_date,
                                 call = rlang::caller_env(),
                                 add_err_msg = "") {
  checkmate::assert_date(dates1)
  checkmate::assert_date(dates2)
  check_data_overlap <- min(dates1) <= max(dates2) &
    max(dates1) <= max_date
  if (!check_data_overlap) {
    cli::cli_abort(
      c(
        "The two vectors of dates do not overlap",
        add_err_msg
      ),
      call = call,
      class = "wwinference_input_data_error"
    )
  }

  invisible()
}
#' Assert that two tibbles of date and time mapping align
#'
#' @param first_data a tibble containing the columns `date` (with IS08601
#' dates) and `t` (integers of time in days)
#' @param second_data a tibble containing the columns `date` (with
#' IS08601 dates) and `t` (integers of time in days)
#' @param arg string to print the name of the element your checking
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#' @param add_err_msg add_err_msg string containing an additional error message,
#' default is the empty string (`""`)
#'
#' @return NULL invisible
assert_equivalent_indexing <- function(first_data,
                                       second_data,
                                       arg = "x",
                                       call = rlang::caller_env(),
                                       add_err_msg = "") {
  full_dates <- first_data |> dplyr::distinct(date, t)
  subset_dates <- second_data |>
    dplyr::distinct(date, t) |>
    dplyr::rename(second_t = t)

  test_df <- full_dates |>
    dplyr::left_join(subset_dates, by = "date") |>
    dplyr::filter(!is.na(second_t))

  check_indexing <- all(test_df$t == test_df$second_t)

  if (!check_indexing) {
    cli::cli_abort(
      c(
        "Date and time indexing on {.arg {arg}} being do not align"
      ),
      call = call,
      class = "wwinference_preprocessing_error"
    )
  }

  invisible()
}
