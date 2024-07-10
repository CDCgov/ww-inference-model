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
#' @param fit_obj a CmdStan object that is the output of fitting the model to
#' the `ww_data` and `count_data`
#' @param date_time_spine A tibble mapping the time index in stan (observed +
#' nowcast + forecast) to real dates
#' @param lab_site_spine A tibble mapping the site-lab index in stan to the
#' corresponding site, lab, and site population
#' @param subpop_spine A tibble mapping the site index in stan to the
#' corresponding subpopulation (either a site or the auxiliary site we add to
#' represent the rest of the population)
#'
#' @return  A tibble containing the full set of posterior draws of the
#' estimated, nowcasted, and forecasted: counts, site-level wastewater
#' concentrations, "global"(e.g. state) R(t) estimate, and the  "local" (site +
#' the one auxiliary subpopulation) R(t) estimates. In the instance where there
#' are observations, the data will be joined to each draw of the predicted
#' observation to facilitate plotting.
#' @export
get_draws_df <- function(ww_data,
                         count_data,
                         fit_obj,
                         date_time_spine,
                         lab_site_spine,
                         subpop_spine) {
  draws <- fit_obj$result$draws()

  count_draws <- draws |>
    tidybayes::spread_draws(pred_hosp[t]) |>
    dplyr::rename(pred_value = pred_hosp) |>
    dplyr::mutate(
      draw = `.draw`,
      name = "pred_counts"
    ) |>
    dplyr::select(name, t, pred_value, draw) |>
    dplyr::left_join(date_time_spine, by = "t") |>
    dplyr::left_join(count_data, by = "date") |>
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
    dplyr::left_join(ww_data, by = c(
      "lab_site_index", "date",
      "lab", "site", "site_pop"
    )) |>
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
    dplyr::left_join(count_data, by = "date") |>
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
