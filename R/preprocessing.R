#' Pre-process wastewater input data, adding needed indices and flagging
#' potential outliers
#' @param ww_data dataframe containing the following columns: site, lab,
#' date, a column for concentration, and lod
#' @param conc_col_name string indicating the name of the column containing
#' the concentration measurements in the wastewater data, default is
#'  `genome_copies_per_ml`
#' @param lod_col_name string indicating the name of the column containing
#' the concentration measurements in the wastewater data, default is
#'  `genome_copies_per_ml`
#' @return a dataframe containing the same columns as ww_data except
#' the `conc_col_name` will be replaced with `genome_copies_per_ml` and
#' the `lod_col_name` will be replaced with `lod` plus the following
#' additional columns needed for the stan model:
#' lab_site_index, site_index, flag_as_ww_outlier, lab_site_name,
#' forecast_date
#' @export
#'
#' @examples
#' ww_data <- tibble::tibble(
#'   date = rep(c("2023-11-01", "2023-11-02"), 2),
#'   site = c(rep(1, 2), rep(2, 2)),
#'   lab = c(1, 1, 1, 1),
#'   conc = c(345.2, 784.1, 401.5, 681.8),
#'   lod = c(20, 20, 15, 15)
#' )

#' ww_data_preprocessed <- preprocess_ww_data(ww_data,
#'                                            conc_col_name = "conc",
#'                                            lod_col_name = "lod"
#'                                            )
preprocess_ww_data <- function(ww_data,
                               conc_col_name = "genome_copies_per_ml",
                               lod_col_name = "lod") {
  # Add some columns
  ww_data_add_cols <- ww_data |>
    dplyr::left_join(
      ww_data |>
        dplyr::distinct(lab, site) |>
        dplyr::mutate(
          lab_site_index = dplyr::row_number()
        ),
      by = c("lab", "site")
    ) |>
    dplyr::left_join(
      ww_data |>
        dplyr::distinct(site) |>
        dplyr::mutate(site_index = dplyr::row_number()),
      by = "site"
    ) |>
    dplyr::rename(
      lod = {{ lod_col_name }},
      genome_copies_per_ml = {{ conc_col_name }}
    ) |>
    dplyr::mutate(
      lab_site_name = glue::glue("Site: {site},  Lab:  {lab}"),
      below_lod = ifelse(genome_copies_per_ml < lod, 1, 0)
    )

  # Get an extra column that identifies the wastewater outliers using the
  # default parameters
  ww_preprocessed <- flag_ww_outliers(ww_data_add_cols,
    conc_col_name = "genome_copies_per_ml"
  )

  return(ww_preprocessed)
}


#' Pre-process hospital admissions data, converting column names to those
#' that [get_stan_data()] expects.
#' @param hosp_data dataframe containing the following columns: date,
#' a count column, and a population size column
#' @param count_col_name name of the column containing the epidemiological
#' indicator, default is `daily_hosp_admits`
#' @param pop_size_col_name name of the column containing the population size
#' of that the counts are coming from, default is `state_pop`
#'
#' @return a dataframe containing the hospital admissions data renamed to
#' have the following columns `date`, `count`, and `total_pop`
#' @export
#'
#' @examples
#' hosp_data <- tibble::tibble(
#'   date = c("2023-11-01", "2023-11-02"),
#'   daily_admits = c(10, 20),
#'   state_pop = c(1e6, 1e6)
#' )
#' hosp_data_preprocessed <- preprocess_hosp_data(
#'   hosp_data,
#'   "daily_admits",
#'   "state_pop"
#' )
preprocess_hosp_data <- function(hosp_data,
                                 count_col_name = "daily_hosp_admits",
                                 pop_size_col_name = "state_pop") {
  hosp_data_preprocessed <- hosp_data |>
    dplyr::rename(
      count = {{ count_col_name }},
      total_pop = {{ pop_size_col_name }}
    )

  return(hosp_data_preprocessed)
}




#' Flag WW outliers
#'
#' @param ww_data dataframe containing the following columns: site, lab,
#' lab_site_index, date, a column for concentration, and below_lod
#' @param conc_col_name string, name of the column containing the concentration
#' measurements in the wastewater data, default is `genome_copies_per_ml`
#' @param rho_threshold float indicating the z-score threshold for "jump"
#' @param log_conc_threshold float indicating the z-score threshold for
#' log concentration
#' @param threshold_n_dps min number of data points above the LOD per lab-site
#'
#' @return ww_w_outliers_flaged dataframe containing all of the columns in
#' ww_data input dataframe plus two additional columns:
#'  `flag_as_ww_outlier` and `exclude`
#' `flag as_ww_outlier` contains a 0 if the datapoint is not an outlier and a 1
#' if it is an outlier. `exclude` tells the model whether or not to exclude that
#' data point, which here is by default set to 0 for all data points (even
#' those flagged as outliers). Excluding the outliers is a second optional
#' step.
#' @export
#'
#' @examples
#' ww_data <- wwinference::ww_data
#' ww_data_preprocessed <- wwinference::preprocess_ww_data(ww_data)
#' ww_data_outliers_flagged <- flag_ww_outliers(ww_data_preprocessed)
flag_ww_outliers <- function(ww_data,
                             conc_col_name = "genome_copies_per_ml",
                             rho_threshold = 2,
                             log_conc_threshold = 3,
                             threshold_n_dps = 1) {
  n_dps <- ww_data |>
    dplyr::filter(below_lod == 0) |>
    dplyr::group_by(lab_site_index) |>
    dplyr::summarise(n_data_points = dplyr::n())

  # Get the ww statistics we need for outlier detection
  ww_stats <- ww_data |>
    dplyr::left_join(n_dps,
      by = "lab_site_index"
    ) |>
    # exclude below LOD from z scoring and remove lab-sites with too
    # few data points
    dplyr::filter(
      below_lod == 0,
      n_data_points > threshold_n_dps
    ) |>
    dplyr::group_by(lab_site_index) |>
    dplyr::arrange(date, "desc") |>
    dplyr::mutate(
      log_conc = log(!!rlang::sym(conc_col_name)),
      prev_log_conc = dplyr::lag(log_conc, 1),
      prev_date = dplyr::lag(date, 1),
      diff_log_conc = log_conc - prev_log_conc,
      diff_time = as.numeric(difftime(date, prev_date)),
      rho = diff_log_conc / diff_time
    ) |>
    dplyr::select(date, lab_site_index, rho) |>
    dplyr::distinct()

  # Combine stats with ww data
  ww_rho <- ww_data |>
    left_join(ww_stats, by = c("lab_site_index", "date"))

  # compute z scores and flag
  ww_z_scored <- ww_rho |>
    dplyr::left_join(
      ww_rho |>
        dplyr::group_by(lab_site_index) |>
        dplyr::summarise(
          mean_rho = mean(rho, na.rm = TRUE),
          std_rho = sd(rho, na.rm = TRUE),
          mean_conc = mean(!!rlang::sym(conc_col_name), na.rm = TRUE),
          std_conc = sd(!!rlang::sym(conc_col_name), na.rm = TRUE)
        ),
      by = "lab_site_index"
    ) |>
    dplyr::group_by(lab_site_index) |>
    mutate(
      z_score_conc = (!!rlang::sym(conc_col_name) - mean_conc) / std_conc,
      z_score_rho = (rho - mean_rho) / std_rho
    ) |>
    dplyr::mutate(
      z_score_rho_t_plus_1 = dplyr::lead(z_score_rho, 1),
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
    dplyr::ungroup() |>
    dplyr::mutate(
      exclude = 0 # by default, we don't exclude anything
    )

  ww_w_outliers_flagged <- ww_z_scored |>
    dplyr::select(
      colnames(ww_data),
      flag_as_ww_outlier,
      exclude
    )

  return(ww_w_outliers_flagged)
}

#' Indicate data that we want to exclude from model fitting
#' @description This function takes in a dataframe which contains an outlier
#' column name specified by the `outlier_col_name`.
#'
#' @param data A dataframe of preprocessed data to be used to fit the o
#' @param outlier_col_name A character string indicating the name of the column
#' containing the outlier indicator, must contain only 0 or 1
#' @param remove_outliers A boolean indicating whether or not to exclude the
#' outliers from the fitting. If TRUE, copy outliers to exclusions, if FALSE,
#' set exclusions to none
#'
#' @return a dataframe with the same columns as in `data` plus an additional
#' `exclude` column containing 0s for the data to be passed to the model
#' and 1s where the data should be excluded
#' @export
#'
#' @examples
#' data <- tibble::tibble(
#'   date = c("2023-10-01", "2023-10-02"),
#'   genome_copies_per_mL = c(300, 3e6),
#'   flag_as_ww_outlier = c(0, 1),
#'   exclude = c(0, 0)
#' )
#' data_w_exclusions <- indicate_ww_exclusions(data,
#'   outlier_col_name = "flag_as_ww_outlier",
#'   remove_outliers = TRUE
#' )
indicate_ww_exclusions <- function(data,
                                   outlier_col_name = "flag_as_ww_outlier",
                                   remove_outliers = TRUE) {
  # Check for the presence of the outlier column name
  if (!outlier_col_name %in% c(colnames(data))) {
    cli::cli_abort(
      "Specified name of the outlier column not present in the data"
    )
  }


  if (isTRUE(remove_outliers)) {
    # Port over the outliers flagged to the exclude column
    data_w_exclusions <- data |>
      dplyr::mutate(
        exclude = ifelse({{ outlier_col_name }} == 1, 1, exclude)
      )
  }
  return(data_w_exclusions)
}
