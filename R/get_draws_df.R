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
#' @param count_data A dataframe of the preprocessed daily count data (e.g.
#' hospital admissions) from the "global" population
#' @param stan_data_list A list containing all the data passed to stan for
#' fitting the model
#' @param fit_obj a CmdStan object that is the output of fitting the model to
#' `x` and `count_data`
#' @param ... additional arguments
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

#' @rdname get_draws
#' @export
get_draws.data.frame <- function(x,
                                 count_data,
                                 stan_data_list,
                                 fit_obj,
                                 ...,
                                 what = "all") {
  # Checking we are getting all
  what_ok <- c(
    "all", "predicted_counts", "predicted_ww", "global_rt", "site_level_rt"
  )

  if (any(!what %in% what_ok)) {
    idx <- which(!what %in% what_ok)
    stop(
      "The following invalid values were passed to `what`: ",
      paste(what[idx], collapse = ", "), ". Valid values incllude: ",
      paste(what_ok, collapse = ", "), "."
    )
  }

  names(what_ok) <- what_ok
  what_ok <- FALSE
  if ("all" %in% what) {
    what_ok[] <- TRUE
  } else {
    what_ok[what] <- TRUE
  }

  draws <- fit_obj$result$draws()

  # Get the necessary mappings needed to join draws to data
  date_time_spine <- tibble::tibble(
    date = seq(
      from = min(count_data$date),
      to = min(count_data$date) + stan_data_list$ot + stan_data_list$ht,
      by = "days"
    )
  ) |>
    dplyr::mutate(t = row_number())

  # Lab-site index to corresponding lab, site, and site population size
  lab_site_spine <- x |>
    dplyr::distinct(.data$site, .data$lab, .data$lab_site_index, .data$site_pop)

  # Site index to corresponding site and subpopulation size
  subpop_spine <- x |>
    dplyr::distinct(.data$site, .data$site_index, .data$site_pop) |>
    dplyr::mutate(site = as.factor(.data$site)) |>
    dplyr::bind_rows(tibble::tibble(
      site = "remainder of pop",
      site_index = max(x$site_index) + 1,
      site_pop = stan_data_list$subpop_size[
        length(unique(stan_data_list$subpop_size))
      ]
    ))

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
      dplyr::select(-"t")
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
      dplyr::left_join(lab_site_spine, by = "lab_site_index") |>
      dplyr::left_join(
        x |>
          dplyr::select(-"t"),
        by = c(
          "lab_site_index", "date",
          "lab", "site", "site_pop"
        )
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(observed_value = .data$log_genome_copies_per_ml) |>
      dplyr::select(-"t")
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
      dplyr::rename("observed_value" = "count") |>
      dplyr::select(-"t")
  } else {
    NULL
  }

  site_level_rt_draws <- if (what_ok["site_level_rt"]) {
    draws |>
      tidybayes::spread_draws(!!str2lang("r_site_t[site_index, t]")) |>
      dplyr::rename("pred_value" = "r_site_t") |>
      dplyr::mutate(
        draw = .data$`.draw`,
        pred_value = .data$pred_value
      ) |>
      dplyr::select("site_index", "t", "pred_value", "draw") |>
      dplyr::left_join(date_time_spine, by = "t") |>
      dplyr::left_join(subpop_spine, by = "site_index") |>
      dplyr::ungroup() |>
      dplyr::select(-"t")
  } else {
    NULL
  }

  return(
    new_wwinference_fit_draws(
      predicted_counts = count_draws,
      predicted_ww = ww_draws,
      global_rt = global_rt_draws,
      site_level_rt = site_level_rt_draws
    )
  )
}

#' @export
print.wwinference_fit_draws <- function(x, ...) {
  cat("Draws from the model featuring the following datasets:\n")

  if (length(x$predicted_counts)) {
    cat(
      sprintf(
        "- predicted_counts with %i observations.\n",
        nrow(x$predicted_counts)
      )
    )
  }
  if (length(x$predicted_ww)) {
    cat(
      sprintf(
        "- predicted_ww with %i observations.\n",
        nrow(x$predicted_ww)
      )
    )
  }
  if (length(x$global_rt)) {
    cat(
      sprintf(
        "- global_rt with %i observations.\n",
        nrow(x$global_rt)
      )
    )
  }
  if (length(x$site_level_rt)) {
    cat(
      sprintf(
        "- site_level_rt with %i observations.\n",
        nrow(x$site_level_rt)
      )
    )
  }

  cat("You can use $ to access the datasets.")

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
    site_level_rt) {
  # Checking colnames: Must match all exactly
  predicted_counts_cnames <- c(
    "date", "draw", "observed_value", "pred_value", "total_pop"
  )
  if (length(predicted_counts)) {
    stopifnot(all(sort(colnames(predicted_counts)) == predicted_counts_cnames))
  }

  predicted_ww_cnames <- c(
    "below_lod", "date", "draw", "exclude", "flag_as_ww_outlier",
    "lab", "lab_site_index", "lab_site_name", "log_genome_copies_per_ml",
    "log_lod", "observed_value", "pred_value", "site", "site_index",
    "site_pop"
  )
  if (length(predicted_ww)) {
    stopifnot(all(sort(colnames(predicted_ww)) == predicted_ww_cnames))
  }

  global_rt_cnames <- c(
    "date", "draw", "observed_value", "pred_value", "total_pop"
  )
  if (length(global_rt)) {
    stopifnot(all(sort(colnames(global_rt)) == global_rt_cnames))
  }

  site_level_rt_cnamss <- c(
    "date", "draw", "pred_value", "site", "site_index", "site_pop"
  )
  if (length(site_level_rt)) {
    stopifnot(all(sort(colnames(site_level_rt)) == site_level_rt_cnamss))
  }

  structure(
    list(
      predicted_counts = predicted_counts,
      predicted_ww = predicted_ww,
      global_rt = global_rt,
      site_level_rt = site_level_rt
    ),
    class = "wwinference_fit_draws"
  )
}
