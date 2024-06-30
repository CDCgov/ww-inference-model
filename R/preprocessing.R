#' Get input wastewater data
#' @param ww_data dataframe containing the following columns: site, lab,
#' date, a column for concentration, and lod
#' @param forecast_date The forecast date for this iteration,
#' formatted as a character string in IS08601 format (YYYY-MM-DD).
#' @param conc_col_name name of the column containing the concentration
#' measurements in the wastewater data, default is `genome_copies_per_ml`
#'
#' @return a dataframe containing the transformed and clean NWSS data
#' at the site and lab label for the forecast date and location specified
#' @export
#'
#' @examples
#' ww_data_preprocessed <- preprocess_ww_data(ww_data, "2023-12-01")
preprocess_ww_data <- function(ww_data,
                               forecast_date,
                               conc_col_name = "genome_copies_per_ml") {
  # Add some columns
  ww_data <- ww_data |>
    dplyr::left_join(
      ww_data |>
        dplyr::distinct(lab, site) |>
        dplyr::mutate(
          lab_site = dplyr::row_number()
        ),
      by = c("lab", "site")
    ) |>
    dplyr::mutate(
      lab_site_name = glue::glue("Site: {site},  Lab:  {lab}"),
      below_lod = ifelse({{ conc_col_name }} < lod, 1, 0)
    )

  # Get an extra column that identifies the wastewater outliers using the
  # default parameters
  ww_w_outliers <- flag_ww_outliers(ww_data) |>
    dplyr::mutate(
      forecast_date = !!forecast_date
    ) |>
    # In case the wastewater data being passed in isn't vintaged, we want to
    # make sure we don't include values that are past the forecast date
    dplyr::filter(
      date < forecast_date
    )


  return(ww_w_outliers)
}



#' Flag WW outliers
#'
#' @param ww_data dataframe containing the following columns: site, lab,
#' lab_site, date, a column for concentration, and below_lod
#' @param conc_col_name string, name of the column containing the concentration
#' measurements in the wastewater data, default is `genome_copies_per_ml`
#' @param rho_threshold float indicating the z-score threshold for "jump"
#' @param log_conc_threshold float indicating the z-score threshold for
#' log concentration
#' @param threshold_n_dps min number of data points above the LOD per lab-site
#'
#' @return ww_w_outliers_flaged dataframe containing all of the columns in
#' ww_data input dataframe plus an additional column `flag_as_ww_outlier`
#' which contains a 0 if the datapoint is not an outlier and a 1 if it is
#' an outlier.
#' @export
#'
#' @examples
#' ww_data_outliers_flagged <- flag_ww_outliers(ww_data)
flag_ww_outliers <- function(ww_data,
                             conc_col_name = "genome_copies_per_ml",
                             rho_threshold = 2,
                             log_conc_threshold = 3,
                             threshold_n_dps = 1) {
  n_dps <- ww_data |>
    dplyr::filter(below_lod == 0) |>
    dplyr::group_by(lab_site) |>
    dplyr::summarise(n_data_points = dplyr::n())

  # Get the ww statistics we need for outlier detection
  ww_stats <- ww_data |>
    dplyr::left_join(n_dps,
      by = "lab_site"
    ) |>
    # exclude below LOD from z scoring and remove lab-sites with too
    # few data points
    dplyr::filter(
      below_lod == 0,
      n_data_points > threshold_n_dps
    ) |>
    dplyr::group_by(lab_site) |>
    dplyr::arrange(date, "desc") |>
    dplyr::mutate(
      log_conc = log(!!sym(conc_col_name)),
      prev_log_conc = lag(log_conc, 1),
      prev_date = lag(date, 1),
      diff_log_conc = log_conc - prev_log_conc,
      diff_time = as.numeric(difftime(date, prev_date)),
      rho = diff_log_conc / diff_time
    ) |>
    dplyr::select(date, lab_site, rho) |>
    dplyr::distinct()

  # Combine stats with ww data
  ww_rho <- ww_data |>
    left_join(ww_stats, by = c("lab_site", "date"))

  # compute z scores and flag
  ww_z_scored <- ww_rho |>
    dplyr::left_join(
      ww_rho |>
        dplyr::group_by(lab_site) |>
        dplyr::summarise(
          mean_rho = mean(rho, na.rm = TRUE),
          std_rho = sd(rho, na.rm = TRUE),
          mean_conc = mean(!!sym(conc_col_name), na.rm = TRUE),
          std_conc = sd(!!sym(conc_col_name), na.rm = TRUE)
        ),
      by = "lab_site"
    ) |>
    dplyr::group_by(lab_site) |>
    mutate(
      z_score_conc = (!!sym(conc_col_name) - mean_conc) / std_conc,
      z_score_rho = (rho - mean_rho) / std_rho
    ) |>
    dplyr::mutate(
      z_score_rho_t_plus_1 = lead(z_score_rho, 1),
      flagged_for_removal_conc = dplyr::case_when(
        abs(z_score_conc) >= log_conc_threshold ~ 1,
        is.na(z_score_conc) ~ 0,
        TRUE ~ 0
      ),
      flagged_for_removal_rho = dplyr::case_when(
        (
          abs(z_score_rho) >= rho_threshold &
            (abs(z_score_rho_t_plus_1) >= rho_threshold) &
            sign(z_score_rho != sign(z_score_rho_t_plus_1))
        ) ~ 1,
        is.na(z_score_rho) ~ NA,
        TRUE ~ 0
      )
    ) |>
    dplyr::mutate(flag_as_ww_outlier = dplyr::case_when(
      flagged_for_removal_rho == 1 ~ 1,
      flagged_for_removal_conc == 1 ~ 1,
      TRUE ~ 0
    )) |>
    dplyr::ungroup()

  ww_w_outliers_flagged <- ww_z_scored |>
    dplyr::select(
      colnames(ww_data),
      flag_as_ww_outlier
    )

  return(ww_w_outliers_flagged)
}
