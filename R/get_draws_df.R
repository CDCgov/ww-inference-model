#' @title Postprocess to generate a draws dataframe
#'
#' @description
#' This function takes in the two input data sources, the CmdStan fit object,
#' and the 3 relevant mappings from stan indices to the real data, in order
#' to generate a dataframe containing the posterior draws of the counts (e.g.
#' hospital admissions), the wastewater concentration values, the "global" R(t),
#' and the "local" R(t) estimates + the critical metadata in the data
#'
#'
#' @param ww_data A dataframe of the preprocessed wastewater concentration data
#' used to fit the model
#' @param count_data A dataframe of the preprocessed daily count data (e.g.
#' hospital admissions) from the "global" population
#' @param stan_args A list containing all the data passed to stan for fitting
#' the model
#' @param fit_obj a CmdStan object that is the output of fitting the model to
#' the `ww_data` and `count_data`
#' @param ww_output an object of the `wwinference_fit` class
#' @return  A tibble containing the full set of posterior draws of the
#' estimated, nowcasted, and forecasted: counts, site-level wastewater
#' concentrations, "global"(e.g. state) R(t) estimate, and the  "local" (site +
#' the one auxiliary subpopulation) R(t) estimates. In the instance where there
#' are observations, the data will be joined to each draw of the predicted
#' observation to facilitate plotting.
#' @export
get_draws_df <- function(ww_output, ...) UseMethod("get_draws_df")

#' @export
#' @rdname get_draws_df
get_draws_df.wwinference_fit <- function(ww_output, ...) {
  get_draws_df.default(
    ww_data = ww_output$input_data$input_ww_data,
    count_data = ww_output$input_data$input_count_data,
    stan_args = ww_output$stan_args,
    fit_obj = ww_output$fit
  )
}

#' @export
#' @rdname get_draws_df
get_draws_df.default <- function(ww_data,
                                 count_data,
                                 stan_args,
                                 fit_obj) {
  draws <- fit_obj$result$draws()

  # Get the necessary mappings needed to join draws to data
  date_time_spine <- tibble::tibble(
    date = seq(
      from = min(count_data$date),
      to = min(count_data$date) + stan_args$ot + stan_args$ht,
      by = "days"
    )
  ) |>
    dplyr::mutate(t = row_number())
  # Lab-site index to corresponding lab, site, and site population size
  lab_site_spine <- ww_data |>
    dplyr::distinct(site, lab, lab_site_index, site_pop)
  # Site index to corresponding site and subpopulation size
  subpop_spine <- ww_data |>
    dplyr::distinct(site, site_index, site_pop) |>
    dplyr::mutate(site = as.factor(site)) |>
    dplyr::bind_rows(tibble::tibble(
      site = "remainder of pop",
      site_index = max(ww_data$site_index) + 1,
      site_pop = stan_args$subpop_size[
        length(unique(stan_args$subpop_size))
      ]
    ))


  count_draws <- draws |>
    tidybayes::spread_draws(pred_hosp[t]) |>
    dplyr::rename(pred_value = pred_hosp) |>
    dplyr::mutate(
      draw = `.draw`,
      name = "pred_counts"
    ) |>
    dplyr::select(name, t, pred_value, draw) |>
    dplyr::left_join(date_time_spine, by = "t") |>
    dplyr::left_join(
      count_data |>
        dplyr::select(-t),
      by = "date"
    ) |>
    dplyr::ungroup() |>
    dplyr::rename(observed_value = count) |>
    dplyr::mutate(
      observation_type = "count",
      type_of_quantity = "global",
      lab_site_index = NA,
      subpop = NA,
      lab = NA,
      site_pop = NA,
      below_lod = NA,
      lod = NA,
      flag_as_ww_outlier = NA,
      exclude = NA
    ) |>
    dplyr::select(-t)

  ww_draws <- draws |>
    tidybayes::spread_draws(pred_ww[lab_site_index, t]) |>
    dplyr::rename(pred_value = pred_ww) |>
    dplyr::mutate(
      draw = `.draw`,
      name = "pred_ww",
      pred_value = exp(pred_value)
    ) |>
    dplyr::select(name, lab_site_index, t, pred_value, draw) |>
    dplyr::left_join(date_time_spine, by = "t") |>
    dplyr::left_join(lab_site_spine, by = "lab_site_index") |>
    dplyr::left_join(
      ww_data |>
        dplyr::select(-t),
      by = c(
        "lab_site_index", "date",
        "lab", "site", "site_pop"
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(observed_value = genome_copies_per_ml) |>
    dplyr::mutate(
      observation_type = "genome copies per mL",
      type_of_quantity = "local",
      total_pop = NA,
      subpop = glue::glue("Site: {site}")
    ) |>
    dplyr::select(colnames(count_draws), -t)

  global_rt_draws <- draws |>
    tidybayes::spread_draws(rt[t]) |>
    dplyr::rename(pred_value = rt) |>
    dplyr::mutate(
      draw = `.draw`,
      name = "global R(t)"
    ) |>
    dplyr::select(name, t, pred_value, draw) |>
    dplyr::left_join(date_time_spine, by = "t") |>
    dplyr::left_join(
      count_data |>
        dplyr::select(-t),
      by = "date"
    ) |>
    dplyr::ungroup() |>
    dplyr::rename(observed_value = count) |>
    dplyr::mutate(
      observed_value = NA,
      observation_type = "latent variable",
      type_of_quantity = "global",
      lab_site_index = NA,
      subpop = NA,
      lab = NA,
      site_pop = NA,
      below_lod = NA,
      lod = NA,
      flag_as_ww_outlier = NA,
      exclude = NA
    ) |>
    dplyr::select(-t)

  site_level_rt_draws <- draws |>
    tidybayes::spread_draws(r_site_t[site_index, t]) |>
    dplyr::rename(pred_value = r_site_t) |>
    dplyr::mutate(
      draw = `.draw`,
      name = "subpop R(t)",
      pred_value = pred_value
    ) |>
    dplyr::select(name, site_index, t, pred_value, draw) |>
    dplyr::left_join(date_time_spine, by = "t") |>
    dplyr::left_join(subpop_spine, by = "site_index") |>
    dplyr::ungroup() |>
    dplyr::mutate(
      observed_value = NA,
      lab_site_index = NA,
      lab = NA,
      below_lod = NA,
      lod = NA,
      flag_as_ww_outlier = NA,
      exclude = NA,
      observation_type = "latent variable",
      type_of_quantity = "local",
      total_pop = NA,
      subpop = ifelse(site != "remainder of pop",
        glue::glue("Site: {site}"), "remainder of pop"
      )
    ) |>
    dplyr::select(colnames(count_draws), -t)

  draws_df <- dplyr::bind_rows(
    count_draws,
    ww_draws,
    global_rt_draws,
    site_level_rt_draws
  )

  return(draws_df)
}
