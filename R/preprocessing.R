#' Get input wastewater data
#' @param ww_data
#' @param conc_col_name name of the column containing the concentration
#' measurements in the wastewater data, default is `genome_copies_per_ml`
#' @param forecast_date The forecast date for this iteration,
#' formatted as a character string in IS08601 format (YYYY-MM-DD).
#' @param calibration_time The duration of the model calibration period
#' (relative to the last hospital admissions data point) in units of
#' model timesteps (typically days).
#' @param last_hosp_data_date A date indicating the date of last reported
#' hospital admission as of the forecast date
#'
#' @return a dataframe containing the transformed and clean NWSS data
#' at the site and lab label for the forecast date and location specified
#' @export
preprocess_ww_data <- function(ww_data,
                               conc_col_name = "genome_copies_per_ml",
                               forecast_date,
                               calibration_time,
                               last_hosp_data_date,
                               ww_data_mapping) {
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

  # Get extra columns that identify wastewater outliers
  ww_w_outliers <- flag_ww_outliers(ww_data) |>
    select(
      date, location, ww, site, lab, lab_wwtp_unique_id,
      ww_pop, below_LOD, lod_sewage, flag_as_ww_outlier
    )
  # If more than one location, than this data isn't being used for fitting
  # And we don't wanto generate these
  if (length(location_i) == 1) {
    site_map <- ww_w_outliers |>
      distinct(site) |>
      mutate(site_index = row_number())
    site_lab_map <- ww_w_outliers |>
      distinct(lab_wwtp_unique_id) |>
      mutate(lab_site_index = row_number())

    ww <- ww_w_outliers |>
      left_join(site_map, by = "site") |>
      left_join(site_lab_map, by = "lab_wwtp_unique_id")
  } else {
    ww <- ww_w_outliers
  }


  return(ww)
}



#' Flag WW outliers
#'
#' @param ww_data data at the lab-site level of WW concentrations
#' @param rho_threshold z-score threshold for "jump"
#' @param log_conc_threshold z-score threshold for log concentration
#' @param threshold_n_dps min number of data points above the LOD per lab-site
#'
#' @return ww_data + columns for outlier flagging
#' @export
#'
#' @examples
flag_ww_outliers <- function(ww_data,
                             rho_threshold = 2,
                             log_conc_threshold = 3,
                             threshold_n_dps = 1) {
  n_dps <- ww_data |>
    dplyr::filter(below_LOD == 0) %>%
    group_by(lab_wwtp_unique_id) %>%
    summarise(n_data_points = n())

  # Get the ww statistics we need for outlier detection
  ww_stats <- ww_data %>%
    left_join(n_dps, by = "lab_wwtp_unique_id") %>%
    # exclude below LOD from z scoring and remove lab-sites with too
    # few data points
    dplyr::filter(below_LOD == 0, n_data_points > threshold_n_dps) %>%
    group_by(lab_wwtp_unique_id) %>%
    arrange(date, "desc") %>%
    mutate(
      log_conc = log(ww),
      prev_log_conc = lag(log_conc, 1),
      prev_date = lag(date, 1),
      diff_log_conc = log_conc - prev_log_conc,
      diff_time = as.numeric(difftime(date, prev_date)),
      rho = diff_log_conc / diff_time
    ) %>%
    select(date, lab_wwtp_unique_id, rho) %>%
    distinct()

  # Combine stats with ww data
  ww_rho <- ww_data %>%
    left_join(ww_stats, by = c("lab_wwtp_unique_id", "date"))

  # compute z scores and flag
  ww_z_scored <- ww_rho %>%
    dplyr::left_join(
      ww_rho %>%
        dplyr::group_by(lab_wwtp_unique_id) %>%
        dplyr::summarise(
          mean_rho = mean(rho, na.rm = TRUE),
          std_rho = sd(rho, na.rm = TRUE),
          mean_conc = mean(ww, na.rm = TRUE),
          std_conc = sd(ww, na.rm = TRUE)
        ),
      by = "lab_wwtp_unique_id"
    ) %>%
    dplyr::group_by(lab_wwtp_unique_id) %>%
    mutate(
      z_score_conc = (ww - mean_conc) / std_conc,
      z_score_rho = (rho - mean_rho) / std_rho
    ) %>%
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
    ) %>%
    dplyr::mutate(flag_as_ww_outlier = dplyr::case_when(
      flagged_for_removal_rho == 1 ~ 1,
      flagged_for_removal_conc == 1 ~ 1,
      TRUE ~ 0
    )) %>%
    dplyr::ungroup()

  return(ww_z_scored)
}
