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
get_draws_df <- function(x, ...) {
  UseMethod("get_draws_df")
}

#' S3 method for extracting posterior draws alongside data for a
#' wwinference_fit object
#'
#' This method overloads the generic get_draws_df function specifically
#' for objects of type 'wwinference_fit'.
#'
#' @rdname get_draws_df
#' @export
get_draws_df.wwinference_fit <- function(x, ...) {
  get_draws_df.data.frame(
    x = x$raw_input_data$input_ww_data,
    count_data = x$raw_input_data$input_count_data,
    stan_data_list = x$stan_data_list,
    fit_obj = x$fit
  )
}

#' @export
#' @rdname get_draws_df
get_draws_df.default <- function(x, ...) {
  stop(
    "No method defined for get_draws_df for object of class(es) ",
    paste(class(x), collapse = ", "),
    ". Use directly on a wwinference_fit object or a",
    "dataframe of wastewater observations.",
    call. = FALSE
  )
}

#' @rdname get_draws_df
#' @export
get_draws_df.data.frame <- function(x,
                                    count_data,
                                    stan_data_list,
                                    fit_obj,
                                    ...) {
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


  count_draws <- draws |>
    tidybayes::spread_draws(!!str2lang("pred_hosp[t]")) |>
    dplyr::rename("pred_value" = "pred_hosp") |>
    dplyr::mutate(
      draw = .data$`.draw`,
      name = "predicted counts"
    ) |>
    dplyr::select("name", "t", "pred_value", "draw") |>
    dplyr::left_join(date_time_spine, by = "t") |>
    dplyr::left_join(
      count_data |>
        dplyr::select(-"t"),
      by = "date"
    ) |>
    dplyr::ungroup() |>
    dplyr::rename("observed_value" = "count") |>
    dplyr::mutate(
      observation_type = "count",
      type_of_quantity = "global",
      lab_site_index = NA,
      subpop = NA,
      lab = NA,
      site_pop = NA,
      below_lod = NA,
      log_lod = NA,
      flag_as_ww_outlier = NA,
      exclude = NA
    ) |>
    dplyr::select(-"t")

  ww_draws <- draws |>
    tidybayes::spread_draws(!!str2lang("pred_ww[lab_site_index, t]")) |>
    dplyr::rename("pred_value" = "pred_ww") |>
    dplyr::mutate(
      draw = .data$`.draw`,
      name = "predicted wastewater",
    ) |>
    dplyr::select("name", "lab_site_index", "t", "pred_value", "draw") |>
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
    dplyr::mutate(
      observation_type = "log genome copies per mL",
      type_of_quantity = "local",
      total_pop = NA,
      subpop = glue::glue("Site: {site}")
    ) |>
    dplyr::select(colnames(count_draws), -"t")

  global_rt_draws <- draws |>
    tidybayes::spread_draws(!!str2lang("rt[t]")) |>
    dplyr::rename("pred_value" = "rt") |>
    dplyr::mutate(
      draw = .data$`.draw`,
      name = "global R(t)"
    ) |>
    dplyr::select("name", "t", "pred_value", "draw") |>
    dplyr::left_join(date_time_spine, by = "t") |>
    dplyr::left_join(
      count_data |>
        dplyr::select(-"t"),
      by = "date"
    ) |>
    dplyr::ungroup() |>
    dplyr::rename("observed_value" = "count") |>
    dplyr::mutate(
      observed_value = NA,
      observation_type = "latent variable",
      type_of_quantity = "global",
      lab_site_index = NA,
      subpop = NA,
      lab = NA,
      site_pop = NA,
      below_lod = NA,
      log_lod = NA,
      flag_as_ww_outlier = NA,
      exclude = NA
    ) |>
    dplyr::select(-"t")

  site_level_rt_draws <- draws |>
    tidybayes::spread_draws(!!str2lang("r_site_t[site_index, t]")) |>
    dplyr::rename("pred_value" = "r_site_t") |>
    dplyr::mutate(
      draw = .data$`.draw`,
      name = "subpopulation R(t)",
      pred_value = .data$pred_value
    ) |>
    dplyr::select("name", "site_index", "t", "pred_value", "draw") |>
    dplyr::left_join(date_time_spine, by = "t") |>
    dplyr::left_join(subpop_spine, by = "site_index") |>
    dplyr::ungroup() |>
    dplyr::mutate(
      observed_value = NA,
      lab_site_index = NA,
      lab = NA,
      below_lod = NA,
      log_lod = NA,
      flag_as_ww_outlier = NA,
      exclude = NA,
      observation_type = "latent variable",
      type_of_quantity = "local",
      total_pop = NA,
      subpop = ifelse(.data$site != "remainder of pop",
        glue::glue("Site: {site}"), "remainder of pop"
      )
    ) |>
    dplyr::select(colnames(count_draws), -"t")

  all_draws_df <- dplyr::bind_rows(
    count_draws,
    ww_draws,
    global_rt_draws,
    site_level_rt_draws
  )


  return(all_draws_df)
}
