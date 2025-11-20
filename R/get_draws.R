#' @title Postprocess to generate a draws dataframe
#'
#' @description
#' This function takes in the two input data sources, the CmdStan fit object,
#' and the 3 relevant mappings from stan indices to the real data, in order
#' to generate a dataframe containing the posterior draws of the counts (e.g.
#' hospital admissions), the wastewater concentration values, the "global" R(t),
#' and the "local" R(t) estimates + the critical metadata in the data.
#' This funtion has a default method that takes the two sets of input data,
#' the last of stan arguments, and the CmdStan fitting object, as well as an S3
#' method for objects of class 'wwinference_fit'
#'
#'
#' @param x Either a dataframe of wastewater observations, or an object of
#' class wwinference_fit
#' @param ... additional arguments
#' @param what Character vector. Specifies the variables to extract from the
#' draws. It could be any from `"all"` `"predicted_counts"`, `"predicted_ww"`,
#' `"global_rt"`, or `"subpop_rt"`. When `what = "all"` (the default),
#' the function will extract all four variables.
#' @return  A tibble containing the full set of posterior draws of the
#' estimated, nowcasted, and forecasted: counts, site-level wastewater
#' concentrations, "global"(e.g. state) R(t) estimate, and the  "local" (site +
#' the one auxiliary subpopulation) R(t) estimates. In the instance where there
#' are observations, the data will be joined to each draw of the predicted
#' observation to facilitate plotting.
#' @export
get_draws <- function(x, ..., what = "all") {
  UseMethod("get_draws")
}

#' @rdname get_draws
#' @details
#' The function `get_draws_df()` has been deprecated in favor of `get_draws()`.
#'
#' @export
get_draws_df <- function(x, ...) {
  .Deprecated("get_draws")
}

#' S3 method for extracting posterior draws alongside data for a
#' wwinference_fit object
#'
#' This method overloads the generic `get_draws` function specifically
#' for objects of type 'wwinference_fit'.
#'
#' @rdname get_draws
#' @export
get_draws.wwinference_fit <- function(x, ..., what = "all") {
  get_draws.data.frame(
    x = x$raw_input_data$input_ww_data,
    count_data = x$raw_input_data$input_count_data,
    date_time_spine = x$raw_input_data$date_time_spine,
    site_subpop_spine = x$raw_input_data$site_subpop_spine,
    lab_site_subpop_spine = x$raw_input_data$lab_site_subpop_spine,
    stan_data_list = x$stan_data_list,
    fit_obj = x$fit,
    what = what
  )
}

#' @export
#' @rdname get_draws
get_draws.default <- function(x, ..., what = "all") {
  stop(
    "No method defined for get_draws for object of class(es) ",
    paste(class(x), collapse = ", "),
    ". Use directly on a wwinference_fit object or a",
    "dataframe of wastewater observations.",
    call. = FALSE
  )
}

#' Vector of valid values for `what` in `get_draws`
#' @noRd
get_draws_what_ok <- c(
  "all", "predicted_counts", "predicted_ww", "global_rt", "subpop_rt"
)

#' @rdname get_draws
#' @param count_data A dataframe of the preprocessed daily count data (e.g.
#' hospital admissions) from the "global" population
#' @param date_time_spine tibble mapping dates to time in days
#' @param site_subpop_spine tibble mapping sites to subpopulations
#' @param lab_site_subpop_spine tibble mapping lab-sites to subpopulations
#' @param stan_data_list A list containing all the data passed to stan for
#' fitting the model
#' @param fit_obj a CmdStan object that is the output of fitting the model to
#' `x` and `count_data`
#' @export
get_draws.data.frame <- function(x,
                                 count_data,
                                 date_time_spine,
                                 site_subpop_spine,
                                 lab_site_subpop_spine,
                                 stan_data_list,
                                 fit_obj,
                                 ...,
                                 what = "all") {
  # Checking we are getting all
  what_ok <- get_draws_what_ok

  if (any(!what %in% what_ok)) {
    idx <- which(!what %in% what_ok)
    stop(
      "The following invalid values were passed to `what`: ",
      paste(what[idx], collapse = ", "), ". Valid values include: ",
      paste(what_ok, collapse = ", "), "."
    )
  }

  what_ok <- logical(length(what_ok))
  names(what_ok) <- get_draws_what_ok
  what_ok[] <- FALSE
  if ("all" %in% what) {
    if (length(what) > 1) {
      warning("Ignoring other values of `what` when `all` is present.")
    }
    what_ok[] <- TRUE
  } else {
    what_ok[what] <- TRUE
  }
  if (stan_data_list$include_ww == 0) {
    if (any(c("predicted_ww", "subpop_rt") %in% what)) {
      cli::cli_abort(c(
        "Predicted wastewater concentrations and subpopulation R(t)s",
        " can not be returned because the model wasn't fit to ",
        " site-level wastewater data"
      ))
    }
    what_ok["predicted_ww"] <- FALSE
    what_ok["subpop_rt"] <- FALSE
    if (what == "all") {
      warning(c(
        "Model wasn't fit to wastewater data. ",
        "Predicted wastewater concentrations and subpopulation R(t)s",
        "\nestimates will not be returned in the ",
        "`wwinference_fit_draws` object"
      ))
    }
  }

  draws <- fit_obj$result$draws()


  count_draws <- if (what_ok["predicted_counts"]) {
    draws |> # predicted_counts
      tidybayes::spread_draws(!!str2lang("pred_hosp[t]")) |>
      dplyr::rename("pred_value" = "pred_hosp") |>
      dplyr::mutate(
        draw = .data$`.draw`,
      ) |>
      dplyr::select("t", "pred_value", "draw") |>
      dplyr::left_join(date_time_spine, by = "t") |>
      dplyr::left_join(
        count_data |>
          dplyr::select(-"t"),
        by = "date"
      ) |>
      dplyr::ungroup() |>
      dplyr::rename("observed_value" = "count") |>
      dplyr::select(
        "date",
        "draw",
        "observed_value",
        "pred_value",
        "total_pop"
      )
  } else {
    NULL
  }


  ww_draws <- if (what_ok["predicted_ww"]) {
    draws |>
      tidybayes::spread_draws(!!str2lang("pred_ww[lab_site_index, t]")) |>
      dplyr::rename("pred_value" = "pred_ww") |>
      dplyr::mutate(
        draw = .data$`.draw`
      ) |>
      dplyr::select("lab_site_index", "t", "pred_value", "draw") |>
      dplyr::left_join(date_time_spine, by = "t") |>
      dplyr::left_join(lab_site_subpop_spine, by = "lab_site_index") |>
      dplyr::left_join(
        x |> dplyr::distinct(
          .data$log_genome_copies_per_ml,
          .data$log_lod,
          .data$date,
          .data$below_lod,
          .data$lab_site_index
        ),
        by = c(
          "lab_site_index", "date"
        )
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        observed_value = .data$log_genome_copies_per_ml,
      ) |>
      dplyr::select(
        "date",
        "lab_site_name",
        "pred_value",
        "draw",
        "observed_value",
        "subpop_name",
        "subpop_pop",
        "site",
        "lab",
        "log_lod",
        "below_lod",
        "lab_site_index"
      )
  } else {
    NULL
  }

  global_rt_draws <- if (what_ok["global_rt"]) {
    draws |>
      tidybayes::spread_draws(!!str2lang("rt[t]")) |>
      dplyr::rename("pred_value" = "rt") |>
      dplyr::mutate(
        draw = .data$`.draw`
      ) |>
      dplyr::select("t", "pred_value", "draw") |>
      dplyr::left_join(date_time_spine, by = "t") |>
      dplyr::left_join(
        count_data |>
          dplyr::select(-"t"),
        by = "date"
      ) |>
      dplyr::ungroup() |>
      dplyr::select(
        "date",
        "pred_value",
        "draw",
        "total_pop"
      )
  } else {
    NULL
  }

  subpop_rt_draws <- if (what_ok["subpop_rt"]) {
    draws |>
      tidybayes::spread_draws(!!str2lang("r_subpop_t[subpop_index, t]")) |>
      dplyr::rename("pred_value" = "r_subpop_t") |>
      dplyr::mutate(
        draw = .data$`.draw`,
        pred_value = .data$pred_value
      ) |>
      dplyr::select("subpop_index", "t", "pred_value", "draw") |>
      dplyr::left_join(date_time_spine, by = "t") |>
      dplyr::left_join(site_subpop_spine, by = "subpop_index") |>
      dplyr::ungroup() |>
      dplyr::select(
        "date",
        "pred_value",
        "draw",
        "subpop_name",
        "subpop_pop",
      )
  } else {
    NULL
  }

  return(
    new_wwinference_fit_draws(
      predicted_counts = count_draws,
      predicted_ww = ww_draws,
      global_rt = global_rt_draws,
      subpop_rt = subpop_rt_draws
    )
  )
}

#' @export
print.wwinference_fit_draws <- function(x, ...) {
  # Computing the draws
  draws <- c(
    ifelse(length(x$predicted_counts) > 0, max(x$predicted_counts$draw), 0),
    ifelse(length(x$predicted_ww) > 0, max(x$predicted_ww$draw), 0),
    ifelse(length(x$global_rt) > 0, max(x$global_rt$draw), 0),
    ifelse(length(x$subpop_rt) > 0, max(x$subpop_rt$draw), 0)
  ) |> max()

  # This calculates the number of time points in each dataframe
  timepoints <- c(
    ifelse(
      length(x$predicted_counts) > 0,
      diff(range(x$predicted_counts$date)) + 1, 0
    ),
    ifelse(
      length(x$predicted_ww) > 0,
      diff(range(x$predicted_ww$date)) + 1, 0
    ),
    ifelse(
      length(x$global_rt) > 0,
      diff(range(x$global_rt$date)) + 1, 0
    ),
    ifelse(
      length(x$subpop_rt) > 0,
      diff(range(x$subpop_rt$date)) + 1, 0
    )
  ) |> max()

  cat(
    sprintf(
      "Draws from the model featuring %i draws across %i days ",
      draws, timepoints
    ),
    "in the following datasets:\n"
  ) # Same draws and timepoints

  if (length(x$predicted_counts)) {
    cat(
      sprintf(
        " - `$predicted_counts` with %i rows\n",
        nrow(x$predicted_counts)
      )
    )
  }

  if (length(x$predicted_ww)) {
    cat(
      sprintf(
        " - `$predicted_ww` with %i rows across %i sites.\n",
        nrow(x$predicted_ww),
        length(unique(x$predicted_ww$lab_site_index))
      )
    )
  }
  if (length(x$global_rt)) {
    cat(
      sprintf(
        " - `$global_rt` with %i rows\n",
        nrow(x$global_rt)
      )
    )
  }
  if (length(x$subpop_rt)) {
    cat(
      sprintf(
        " - `$subpop_rt` with %i rows across %i subpopulations\n",
        nrow(x$subpop_rt),
        length(unique(x$subpop_rt$subpop_name))
      )
    )
  }

  cat("You can use $ to access the datasets.\n")

  invisible(x)
}

#' Constructor for the new_wwinference_fit_draws
#'
#' Constructor running some checks on the contents of the data.
#'
#' @param predicted_counts Predicted counts
#' @param predicted_ww Predicted ww concentration
#' @param global_rt Global Rt()
#' @param site_level_r Site-level Rt()s
#' @noRd
new_wwinference_fit_draws <- function(
  predicted_counts,
  predicted_ww,
  global_rt,
  subpop_rt
) {
  # Checking colnames: Must match all exactly
  predicted_counts_colnames <- c(
    "date", "pred_value", "observed_value", "draw", "total_pop"
  )
  if (length(predicted_counts)) {
    checkmate::assert_names(
      colnames(predicted_counts),
      permutation.of = predicted_counts_colnames
    )
  }

  predicted_ww_colnames <- c(
    "below_lod",
    "date",
    "draw",
    "lab",
    "lab_site_name",
    "log_lod",
    "observed_value",
    "pred_value",
    "site",
    "subpop_pop",
    "subpop_name",
    "lab_site_index"
  )
  if (length(predicted_ww)) {
    checkmate::assert_names(
      colnames(predicted_ww),
      permutation.of = predicted_ww_colnames
    )
  }

  global_rt_colnames <- c(
    "date", "draw", "pred_value", "total_pop"
  )
  if (length(global_rt)) {
    checkmate::assert_names(
      colnames(global_rt),
      permutation.of = global_rt_colnames
    )
  }

  subpop_rt_colnames <- c(
    "date",
    "draw",
    "pred_value",
    "subpop_pop",
    "subpop_name"
  )
  if (length(subpop_rt)) {
    checkmate::assert_names(
      colnames(subpop_rt),
      permutation.of = subpop_rt_colnames
    )
  }

  structure(
    list(
      predicted_counts = predicted_counts,
      predicted_ww = predicted_ww,
      global_rt = global_rt,
      subpop_rt = subpop_rt
    ),
    class = "wwinference_fit_draws"
  )
}

#' @export
#' @rdname get_draws
#' @param x An object of class `get_draws`.
#' @param y Ignored in the the case of `plot`.
#' @details
#' The plot method for `wwinference_fit_draws` is a wrapper of
#' `get_plot_forecasted_counts`, `get_plot_ww_conc`, `get_plot_global_rt`,
#' and `get_plot_subpop_rt`. Depending on the value of `what`, the function
#' will call the appropriate method.
#'
plot.wwinference_fit_draws <- function(x, y = NULL, what, ...) {
  if (length(what) != 1L) {
    stop(
      "The value provided to `what` must be a length one character vector. ",
      "Currently, it is of length ", length(what), "."
    )
  }

  which_what_are_ok <- setdiff(get_draws_what_ok, "all")

  if (!what %in% which_what_are_ok) {
    stop(
      sprintf(
        paste0(
          "The value provided to what (%s) is invalid. ",
          "Valid values include \"%s\"."
        ),
        paste(what, collapse = ", "),
        paste(which_what_are_ok, collapse = "\", \"")
      )
    )
  }

  if (what == "predicted_counts") {
    get_plot_forecasted_counts(
      draws = x$predicted_counts,
      ...
    )
  } else if (what == "predicted_ww") {
    get_plot_ww_conc(
      x$predicted_ww,
      ...
    )
  } else if (what == "global_rt") {
    get_plot_global_rt(
      x$global_rt,
      ...
    )
  } else if (what == "subpop_rt") {
    get_plot_subpop_rt(
      x$subpop_rt,
      ...
    )
  }
}
