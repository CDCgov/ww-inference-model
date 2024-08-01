#' Check that all dates in dataframe passed in are before a specified date
#'
#' @description
#' This function is specifically meant to ensure that the data in the `df`
#' specified does nto contain data time-stamped beyond the `max_date`. The
#' intended use-case for this is to ensure that one doesn't accidentally
#' pass in data that extends beyond the forecast date, as ideally the user
#' is providing vintaged "as of" datasets or at the very least is filtering the
#' data so that they are not including in their inference data that was
#' made available after the forecast was made.
#'
#'
#' @param df dataframe with `date` column
#' @param max_date string indicating the maximum date in ISO8601 convention
#' e.g. YYYY-MM-DD
#' @param call Calling environment to be passed to the type checker
#'
#' @return NULL, invisibly
check_date_logic <- function(df, max_date, call = rlang::caller_env()) {
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
      c("{.arg {arg}} has negative elements, check to ensure that values",
        "have not been log transformed",
        "!" = "All elements must be 0 or greater",
        "i" = "Elements {.val {which(!is_non_neg)}} are negative"
      ),
      class = "wwinference_input_data_error",
      call = call
    )
  }
  invisible()
}

#' Check that the argument is of type date
#'
#' @param x Object with type checking, here this will be a date
#' @param arg Name of the argument supplying the object
#' @param call Calling environment to be passed to [cli::cli_abort()]
#'
#' @return NULL, invisibly
check_date <- function(x, arg = "x", call = rlang::caller_env()) {
  if ((!rlang::is_integerish(x)) || (!inherits(x, "Date"))) {
    throw_type_error(
      object = x,
      arg_name = arg,
      expected_type = "Date",
      call = call
    )
  }
  invisible()
}

#' Check that the argument is a vector
#'
#' @param x Object with type checking, here this will be a vector
#' @param arg Name of the argument supplying the object
#' @param call Calling environment to be passed to [cli::cli_abort()]
#'
#' @return NULL, invisibly
check_vector <- function(x, arg = "x", call = rlang::caller_env()) {
  # We only want vectors not lists
  if (!rlang::is_bare_vector(x) || inherits(x, "list")) {
    throw_type_error(
      object = x,
      arg_name = arg,
      expected_type = "vector",
      call = call
    )
  }
  invisible()
}


#' Check that the arguments are either integers or characters
#'
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
check_int_or_char <- function(x, arg = "x", call = rlang::caller_env()) {
  # We only want vectors not lists
  if (!(rlang::is_integerish(x) || rlang::is_character(x))) {
    throw_type_error(
      object = x,
      arg_name = arg,
      expected_type = "integer or character",
      call = call
    )
  }
  invisible()
}


#' Check that the arguments are integers
#'
#'
#' @param x Object with type checking
#' @param arg Name of the argument supplying the object
#' @param call Calling environment to be passed to [cli::cli_abort()]
#'
#' @return NULL, invisibly
check_int <- function(x, arg = "x", call = rlang::caller_env()) {
  # We only want vectors not lists
  if (!rlang::is_integerish(x)) {
    throw_type_error(
      object = x,
      arg_name = arg,
      expected_type = "integer",
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
check_required_ww_inputs <- function(ww_data,
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

  # This tells you whats missing
  checkmate::check_names(column_names, must.include = expected_col_names)

  # This tells you from where it is missing
  if (!checkmate::check_names(
    column_names,
    must.include = expected_col_names
  )) {
    cli::cli_abort(
      "Required columns are missing from the input wastewater data",
      class = "wwinference_input_data_error",
      call = call
    )
  }

  invisible()
}

#' Check that the input hosp data contains all the required column names
#'
#' @description
#' This function is intended to be used to check that the hosp data that
#' gets passed into [preprocess_hosp_data()] contains the required columns. If
#' it does not, we want to tell the user which columns are missing. This will
#' not however, ensure that the elements of the column are of the right type,
#' or check that the values of them make sense.
#'
#'
#' @param hosp_data tibble containing the input count data
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
check_required_hosp_inputs <- function(hosp_data,
                                       count_col_name,
                                       pop_size_col_name,
                                       add_req_col_names = c("date"),
                                       call = rlang::caller_env()) {
  column_names <- colnames(hosp_data)
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
  checkmate::check_names(column_names, must.include = expected_col_names)

  # This tells you from where it is missing
  if (!checkmate::check_names(
    column_names,
    must.include = expected_col_names
  )) {
    cli::cli_abort(
      "Required columns are missing from the input count data",
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
  invisible()
}

#' Check that the vector of population sizes for the global catchment area
#' has only a single value
#'
#' @description
#' This function  checks that "global" population sizes in the data passed in as
#' count data contain only a single value of the catchment areas population
#' size. Multiple values might indicate either that there are multiple count
#' data streams being passed in (which is not currently supported) or that
#' there is a time varying population size (which is also not currently
#' supported). This function is specific to the
#' current version of the model, and will
#' be deprecated for a more general `check_for_single_value()` once
#' additional model functionality has been added.
#'
#'
#' @param x the vector to check
#' @param arg the name of the vector to check
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#'
#' @return NULL, invisibly
check_global_pop <- function(x, arg = "x", call = rlang::caller_env()) {
  unique_pops <- unique(x)

  if (length(unique_pops) > 1) {
    cli::cli_abort(
      c("{.arg {arg}} has more than one global population size",
        "i" = c(
          "Multiple/time-varying count catchment area populations",
          "are not currently supported. Check that data is from a single",
          "location, and if so, consider replacing with an average",
          "population size over the inference period"
        )
      ),
      call = call,
      class = "wwinference_input_data_error"
    )
  }
  invisible()
}

#' Check that there are no repeated elements in the vector of dates
#' corresponding to count observations
#'
#' @description
#' This function  checks that the dates in the data passed in as count data are
#' not repeated, which might indicate multiple count data stream. This
#' functionality is not currently supported. This function is specific to the
#' current version of the model, and will
#' be deprecated for a more general `check_for_repeat_elements()` once
#' additional model functionality has been added.
#'
#'
#' @param x the vector to check
#' @param arg the name of the vector to check
#' @param call Calling environment to be passed to [cli::cli_abort()] for
#' traceback.
#'
#' @return NULL, invisibly
check_for_repeat_dates <- function(x, arg = "x", call = rlang::caller_env()) {
  duplicates <- duplicated(x)

  if (sum(duplicates) > 0) {
    cli::cli_abort(
      c("{.arg {arg}} has more than one date",
        "i" = c(
          "Multiple count are not currently supported which might be the",
          "reason there are multiple dates being passed in. ",
          "Check that data is from a single location, and if so, provide",
          "a single count data stream for the population in that location"
        ),
        "!" = "Duplicate element(s) index: {.val {which(duplicates)}}"
      ),
      call = call,
      class = "wwinference_input_data_error"
    )
  }
  invisible()
}
