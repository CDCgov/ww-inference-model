#' Get the input count data to pass directly to stan
#'
#' @param preprocessed_count_data a dataframe with the input count data, must
#' have the following columns: date, count, total_pop
#' @param calibration_time integer indicating the max duration in days that
#' the model is calibrated to the count data for
#' @return datatframe of count data passed to stan
#' @export
get_input_count_data_for_stan <- function(preprocessed_count_data,
                                          calibration_time) {
  # Get the last date that there were observations of the epidemiological
  # indicator (aka cases or hospital admissions counts)
  last_count_data_date <- max(preprocessed_count_data$date, na.rm = TRUE)

  input_count_data_filtered <- preprocessed_count_data |>
    dplyr::filter(
      .data$date > !!last_count_data_date - lubridate::days(!!calibration_time)
    )

  count_data <- add_time_indexing(input_count_data_filtered)


  return(count_data)
}

#' Get the input ww data passed directly to stan
#'
#' @param preprocessed_ww_data a dataframe with the input wastewater data with
#' no gaps, must have the following columns: date, site, lab,
#' genome_copies_per_ml, site_pop, below_lod, and if removing outliers,
#' flag_as_ww_outlier
#' @param first_count_data_date The earliest day with an observation in the '
#' count dataset, in ISO8601 format (YYYY-MM-DD)
#' @param last_count_data_date the last date that a count observation is
#' presen, in ISO8601 format (YYYY-MM-DD)
#' @param calibration_time integer indicating the max duration in days that
#' the model is calibrated to the count data for
#' @return dataframe of the ww data passed to stan
#' @export
get_input_ww_data_for_stan <- function(preprocessed_ww_data,
                                       first_count_data_date,
                                       last_count_data_date,
                                       calibration_time) {
  # Test to see if ww_data_present
  ww_data_present <- !is.null(preprocessed_ww_data)
  if (ww_data_present == FALSE) {
    message("No wastewater data present")
    ww_data <- NULL
  } else {
    if (all(sum(preprocessed_ww_data$flag_as_ww_outlier) > sum(
      preprocessed_ww_data$exclude
    ))) {
      cli::cli_warn(
        c(
          "Wastewater data being passed to the model has outliers flagged,",
          "but not all have been indicated for exclusion from model fit"
        )
      )
    }

    # Test for presence of needed column names
    assert_req_ww_cols_present(preprocessed_ww_data,
      conc_col_name = "log_genome_copies_per_ml",
      lod_col_name = "log_lod"
    )

    # Filter out wastewater outliers, and remove extra wastewater
    # data. Arrange data for indexing. This is what will be returned.
    ww_data <- preprocessed_ww_data |>
      dplyr::filter(
        .data$exclude != 1,
        .data$date > !!last_count_data_date -
          lubridate::days(!!calibration_time)
      ) |>
      dplyr::arrange(.data$date, .data$lab_site_index)
  }

  return(ww_data)
}

#' Get date time spine to map to model output
#'
#' @param forecast_date a character string in ISO8601 format (YYYY-MM-DD)
#' indicating the date that the forecast is to be made.
#' @param input_count_data a dataframe of the count data to be passed
#' directly to stan, , must have the following columns: date, count, total_pop
#' @param last_count_data_date string indicating the date of the last observed
#' count data point in 1SO8601 format (YYYY-MM-DD)
#' @param calibration_time integer indicating the number of days to calibrate
#' the model for, default is `90`
#' @param forecast_horizon integer indicating the number of days, including the
#' forecast date, to produce forecasts for, default is `28`
#'
#'
#' @return a tibble containing an integer for time mapped to the corresponding
#' date, for the entire calibration and forecast period
#' @export
#'
get_date_time_spine <- function(forecast_date,
                                input_count_data,
                                last_count_data_date,
                                calibration_time,
                                forecast_horizon) {
  nowcast_time <- as.integer(
    lubridate::ymd(forecast_date) - last_count_data_date
  )
  date_time_spine <- tibble::tibble(
    date = seq(
      from = min(input_count_data$date),
      to = min(input_count_data$date) +
        calibration_time +
        nowcast_time +
        forecast_horizon,
      by = "days"
    )
  ) |>
    dplyr::mutate(t = row_number())
  return(date_time_spine)
}

#' Get mapping from lab-site to site
#'
#' @param input_ww_data a dataframe of the wastewater data to be passed
#' directly to stan, must have the following columns: date, site, lab,
#' genome_copies_per_ml, site_pop, below_lod, and exclude
#'
#' @return a dataframe mapping the unique combinations of sites and labs
#' to their indices in the model and the population of the site in that
#' observation unit (lab_site)
#' @export
#'
get_lab_site_site_spine <- function(input_ww_data) {
  ww_data_present <- !is.null(input_ww_data)

  if (ww_data_present) {
    lab_site_site_spine <-
      input_ww_data |>
      dplyr::select(
        "lab_site_index", "site_index",
        "site", "lab", "site_pop"
      ) |>
      dplyr::arrange(.data$lab_site_index) |>
      dplyr::distinct() |>
      dplyr::mutate(
        "lab_site_name" = glue::glue(
          "Site: {site}, Lab: {lab}"
        )
      )
  } else {
    lab_site_site_spine <- tibble::tibble()
  }


  return(lab_site_site_spine)
}

#' Get site to subpopulation map
#'
#' @param input_ww_data a dataframe of the wastewater data to be passed
#' directly to stan, must have the following columns: date, site, lab,
#' genome_copies_per_ml, site_pop, below_lod, and exclude
#' @param input_count_data a dataframe of the count data to be passed
#' directly to stan, , must have the following columns: date, count, total_pop
#'
#' @return a dataframe mapping the sites to the corresponding subpopulation and
#' subpopulation index, plus the population in each subpopulation. Imposes
#' the logic to add a subpopulation if the total population is greater than
#' the sum of the site populations in the input wastewater data
#' @export
#'
get_site_subpop_spine <- function(input_ww_data,
                                  input_count_data) {
  ww_data_present <- !is.null(input_ww_data)

  total_pop <- input_count_data |>
    dplyr::distinct(.data$total_pop) |>
    dplyr::pull()

  if (ww_data_present) {
    add_auxiliary_subpop <- ifelse(
      total_pop > sum(unique(input_ww_data$site_pop)),
      TRUE,
      FALSE
    )
    site_indices <- input_ww_data |>
      dplyr::select("site_index", "site", "site_pop") |>
      dplyr::distinct() |>
      dplyr::arrange(.data$site_index)

    if (add_auxiliary_subpop) {
      aux_subpop <- tibble::tibble(
        "site_index" = NA,
        "site" = NA,
        "site_pop" = total_pop - sum(site_indices$site_pop)
      )
    } else {
      aux_subpop <- tibble::tibble()
    }

    site_subpop_spine <- aux_subpop |>
      dplyr::bind_rows(site_indices) |>
      dplyr::mutate(
        subpop_index = dplyr::row_number()
      ) |>
      dplyr::mutate(
        subpop_name = ifelse(!is.na(.data$site),
          glue::glue("Site: {site}"),
          "remainder of population"
        )
      ) |>
      dplyr::rename(
        "subpop_pop" = "site_pop"
      )
  } else {
    site_subpop_spine <- tibble::tibble(
      "site_index" = NA,
      "site" = NA,
      "subpop_pop" = total_pop,
      "subpop_index" = 1,
      "subpop_name" = "total population"
    )
  }

  return(site_subpop_spine)
}

#' Get lab-site subpopulation spine
#'
#' @param lab_site_site_spine tibble mapping lab-sites to sites
#' @param site_subpop_spine tibble mapping sites to subpopulations
#'
#' @return a tibble mapping lab-sites to subpopulations
#' @export
#'
get_lab_site_subpop_spine <- function(lab_site_site_spine,
                                      site_subpop_spine) {
  ww_data_present <- !nrow(lab_site_site_spine) == 0
  # Get lab_site to subpop spine
  if (ww_data_present) {
    lab_site_subpop_spine <- lab_site_site_spine |>
      dplyr::left_join(site_subpop_spine, by = c("site_index", "site"))
  } else {
    lab_site_subpop_spine <- tibble::tibble(
      subpop_index = numeric()
    )
  }

  return(lab_site_subpop_spine)
}


#' Get stan data for ww + hosp model
#'

#' @param input_count_data tibble with the input count data needed for stan
#' @param input_ww_data tibble with the input wastewater data and indices
#' needed for stan
#' @param date_time_spine tibble mapping dates to time in days
#' @param lab_site_site_spine tibble mapping lab-sites to sites
#' @param site_subpop_spine tibble mapping sites to subpopulations
#' @param lab_site_subpop_spine tibble mapping lab-sites to subpopulations
#' @param last_count_data_date string indicating the date of the last data
#' point in the count dataset in ISO8601 convention e.g. YYYY-MM-DD
#' @param first_count_data_date string indicating the date of the first data
#' point in the count dataset in ISO8601 convention e.g. YYYY-MM-DD
#' @param forecast_date string indicating the forecast date in ISO8601
#'  convention e.g. YYYY-MM-DD
#' @param forecast_horizon integer indicating the number of days to make a
#' forecast for
#' @param calibration_time integer indicating the max duration in days that
#' the model is calibrated to the count data for
#' @param generation_interval a vector with a zero-truncated normalized pmf of
#' the generation interval
#' @param inf_to_count_delay a vector with a normalized pmf of the delay from
#'  infection to counts
#' @param infection_feedback_pmf a vector with a normalized pmf dictating the
#' delay of infection feedback
#' @param params a list mapping parameter names to their values
#' @param include_ww integer either 1 or 0 indicating whether to fit the
#' wastewater data or not, default is 1
#' @param compute_likelihood indicator variable telling stan whether or not to
#' compute the likelihood, default = `1`
#'
#' @return `stan_args`: named variables to pass to stan
#' @export
#'
#' @examples
#' ww_data <- tibble::tibble(
#'   date = rep(seq(
#'     from = lubridate::ymd("2023-08-01"),
#'     to = lubridate::ymd("2023-11-01"),
#'     by = "weeks"
#'   ), 2),
#'   site = c(rep(1, 14), rep(2, 14)),
#'   lab = c(rep(1, 28)),
#'   conc = abs(rnorm(28, mean = 500, sd = 50)),
#'   lod = c(rep(20, 14), rep(15, 14)),
#'   site_pop = c(rep(2e5, 14), rep(4e5, 14))
#' )
#'
#' ww_data_preprocessed <- preprocess_ww_data(ww_data,
#'   conc_col_name = "conc",
#'   lod_col_name = "lod"
#' )
#' input_ww_data <- indicate_ww_exclusions(ww_data_preprocessed)
#'
#' hosp_data <- tibble::tibble(
#'   date = seq(
#'     from = lubridate::ymd("2023-07-01"),
#'     to = lubridate::ymd("2023-10-30"),
#'     by = "days"
#'   ),
#'   daily_admits = sample(5:70, 122, replace = TRUE),
#'   state_pop = rep(1e6, 122)
#' )
#'
#' input_count_data <- preprocess_count_data(
#'   hosp_data,
#'   "daily_admits",
#'   "state_pop"
#' )
#'
#' generation_interval <- to_simplex(c(0.01, 0.2, 0.3, 0.2, 0.1, 0.1, 0.01))
#' inf_to_count_delay <- to_simplex(c(
#'   rep(0.01, 12), rep(0.2, 4),
#'   rep(0.01, 10)
#' ))
#' infection_feedback_pmf <- generation_interval
#'
#' params <- get_params(
#'   system.file("extdata", "example_params.toml",
#'     package = "wwinference"
#'   )
#' )
#' forecast_date <- "2023-11-06"
#' calibration_time <- 90
#' forecast_horizon <- 28
#' include_ww <- 1
#' input_count_data_for_stan <- get_input_count_data_for_stan(
#'   input_count_data,
#'   calibration_time
#' )
#' last_count_data_date <- max(input_count_data_for_stan$date, na.rm = TRUE)
#' first_count_data_date <- min(input_count_data_for_stan$date, na.rm = TRUE)
#' input_ww_data_for_stan <- get_input_ww_data_for_stan(
#'   input_ww_data,
#'   first_count_data_date,
#'   last_count_data_date,
#'   calibration_time
#' )
#' date_time_spine <- get_date_time_spine(
#'   forecast_date = forecast_date,
#'   input_count_data = input_count_data_for_stan,
#'   last_count_data_date = last_count_data_date,
#'   forecast_horizon = forecast_horizon,
#'   calibration_time = calibration_time
#' )
#' lab_site_site_spine <- get_lab_site_site_spine(
#'   input_ww_data = input_ww_data_for_stan
#' )
#' site_subpop_spine <- get_site_subpop_spine(
#'   input_ww_data = input_ww_data_for_stan,
#'   input_count_data = input_count_data_for_stan
#' )
#' lab_site_subpop_spine <- get_lab_site_subpop_spine(
#'   lab_site_site_spine = lab_site_site_spine,
#'   site_subpop_spine
#' )
#'
#'
#' stan_data_list <- get_stan_data(
#'   input_count_data_for_stan,
#'   input_ww_data_for_stan,
#'   date_time_spine,
#'   lab_site_site_spine,
#'   site_subpop_spine,
#'   lab_site_subpop_spine,
#'   last_count_data_date,
#'   first_count_data_date,
#'   forecast_date,
#'   forecast_horizon,
#'   calibration_time,
#'   generation_interval,
#'   inf_to_count_delay,
#'   infection_feedback_pmf,
#'   params,
#'   include_ww
#' )
get_stan_data <- function(input_count_data,
                          input_ww_data,
                          date_time_spine,
                          lab_site_site_spine,
                          site_subpop_spine,
                          lab_site_subpop_spine,
                          last_count_data_date,
                          first_count_data_date,
                          forecast_date,
                          forecast_horizon,
                          calibration_time,
                          generation_interval,
                          inf_to_count_delay,
                          infection_feedback_pmf,
                          params,
                          include_ww,
                          compute_likelihood = 1) {
  # Validate input pmfs----------------------------------------------------
  validate_pmf(generation_interval,
    calibration_time,
    input_count_data,
    arg = "generation interval"
  )
  validate_pmf(infection_feedback_pmf,
    calibration_time,
    input_count_data,
    arg = "infection feedback pmf"
  )
  validate_pmf(inf_to_count_delay,
    calibration_time,
    input_count_data,
    arg = "infection to count delay"
  )

  # Check that count data doesn't extend beyond forecast date
  assert_no_dates_after_max(
    date_vector = input_count_data$date,
    max_date = forecast_date,
    arg_dates = "wastewater data",
    arg_max_date = "forecast date"
  )

  # if both datasets are used, validate that that they are
  # compatible and consistent with each other
  if (include_ww == 1) {
    validate_data_jointly(
      input_count_data = input_count_data,
      input_ww_data = input_ww_data,
      date_time_spine = date_time_spine,
      lab_site_site_spine = lab_site_site_spine,
      site_subpop_spine = site_subpop_spine,
      lab_site_subpop_spine = lab_site_subpop_spine,
      calibration_time = calibration_time,
      forecast_date = forecast_date
    )
    # Check that ww data doesn't extend beyond forecast date
    assert_no_dates_after_max(
      date_vector = input_ww_data$date,
      max_date = forecast_date,
      arg_dates = "wastewater data",
      arg_max_date = "forecast date"
    )
  }

  # Define some global variables from the input data-----------------------
  # Get the total pop, coming from the larger population generating the
  # count data
  pop <- input_count_data |>
    dplyr::distinct(.data$total_pop) |>
    dplyr::pull()

  assert_single_value(pop,
    arg = "global population",
    add_err_msg = c(
      "More than one global population size",
      "passed to stan"
    )
  )

  # Get wastewater inputs-------------------------------------------------

  # Returns a list with the numbers of elements needed for the stan model
  ww_data_sizes <- get_ww_data_sizes(
    input_ww_data
  )

  ww_vals <- get_ww_indices_and_values(
    input_ww_data = input_ww_data,
    date_time_spine = date_time_spine,
    lab_site_site_spine = lab_site_site_spine,
    site_subpop_spine = site_subpop_spine,
    lab_site_subpop_spine = lab_site_subpop_spine
  )

  stopifnot(
    "Wastewater sampled times not equal to length of input ww data" =
      length(ww_vals$ww_sampled_times) == ww_data_sizes$owt
  )

  message(
    "Prop of population size covered by wastewater: ",
    sum(unique(input_ww_data$site_pop)) / pop
  )


  # Get count data inputs-----------------------------------------------
  count_data_sizes <- get_count_data_sizes(
    input_count_data = input_count_data,
    forecast_date = forecast_date,
    forecast_horizon = forecast_horizon,
    calibration_time = calibration_time,
    last_count_data_date = last_count_data_date,
    uot = params$uot
  )
  count_indices <- get_count_indices(input_count_data)
  count_values <- get_count_values(
    input_count_data,
    ot = count_data_sizes$ot,
    ht = count_data_sizes$ht,
    count_col_name = "count"
  )

  message(
    "Removed ", nrow(input_ww_data) - ww_data_sizes$owt,
    " outliers from WW data"
  )


  # matrix to transform P(count|I) from weekly to daily
  ind_m <- get_ind_m(
    count_data_sizes$ot + count_data_sizes$ht,
    count_data_sizes$n_weeks
  )
  # matrix to transform p_hosp RW from weekly to daily
  p_hosp_m <- get_ind_m(
    params$uot + count_data_sizes$ot + count_data_sizes$ht,
    count_data_sizes$tot_weeks
  )

  # Estimate of number of initial infections
  i_first_obs_est <- (
    mean(count_values$count[1:7], na.rm = TRUE) /
      params$p_hosp_mean
  )

  # package up parameters for stan data object
  viral_shedding_pars <- c(
    params$t_peak_mean, params$t_peak_sd, params$viral_peak_mean,
    params$viral_peak_sd, params$duration_shedding_mean,
    params$duration_shedding_sd
  )

  inf_to_count_delay_max <- length(inf_to_count_delay)

  stan_data_list <- list(
    gt_max = min(length(generation_interval), params$gt_max),
    hosp_delay_max = inf_to_count_delay_max,
    inf_to_hosp = inf_to_count_delay,
    mwpd = params$ml_of_ww_per_person_day,
    ot = count_data_sizes$ot,
    n_subpops = length(ww_vals$subpop_pops),
    n_ww_sites = ww_data_sizes$n_ww_sites,
    n_ww_lab_sites = ww_data_sizes$n_ww_lab_sites,
    owt = ww_data_sizes$owt,
    oht = count_data_sizes$oht,
    n_censored = ww_data_sizes$n_censored,
    n_uncensored = ww_data_sizes$n_uncensored,
    uot = params$uot,
    ht = count_data_sizes$ht,
    n_weeks = count_data_sizes$n_weeks,
    ind_m = ind_m,
    tot_weeks = count_data_sizes$tot_weeks,
    p_hosp_m = p_hosp_m,
    generation_interval = generation_interval,
    ts = 1:params$gt_max,
    state_pop = pop,
    subpop_size = ww_vals$subpop_pops,
    norm_pop = sum(site_subpop_spine$subpop_pop),
    ww_sampled_times = ww_vals$ww_sampled_times,
    hosp_times = count_indices$count_times,
    ww_sampled_subpops = ww_vals$ww_sampled_subpops,
    lab_site_to_subpop_map = lab_site_subpop_spine$subpop_index,
    ww_sampled_lab_sites = ww_vals$ww_sampled_lab_sites,
    ww_log_lod = ww_vals$ww_lod,
    ww_censored = ww_vals$ww_censored,
    ww_uncensored = ww_vals$ww_uncensored,
    hosp = count_values$counts,
    day_of_week = count_values$day_of_week,
    log_conc = ww_vals$log_conc,
    compute_likelihood = compute_likelihood,
    include_ww = include_ww,
    include_hosp = 1,
    if_l = length(infection_feedback_pmf),
    infection_feedback_pmf = infection_feedback_pmf,
    # All the priors!
    viral_shedding_pars = viral_shedding_pars, # tpeak, viral peak, dur_shed
    autoreg_rt_a = params$autoreg_rt_a,
    autoreg_rt_b = params$autoreg_rt_b,
    autoreg_rt_subpop_a = params$autoreg_rt_subpop_a,
    autoreg_rt_subpop_b = params$autoreg_rt_subpop_b,
    autoreg_p_hosp_a = params$autoreg_p_hosp_a,
    autoreg_p_hosp_b = params$autoreg_p_hosp_b,
    inv_sqrt_phi_prior_mean = params$inv_sqrt_phi_prior_mean,
    inv_sqrt_phi_prior_sd = params$inv_sqrt_phi_prior_sd,
    r_prior_mean = params$r_prior_mean,
    r_prior_sd = params$r_prior_sd,
    log10_g_prior_mean = params$log10_g_prior_mean,
    log10_g_prior_sd = params$log10_g_prior_sd,
    i_first_obs_over_n_prior_a = 1 +
      params$i_first_obs_certainty *
        (i_first_obs_est / pop),
    i_first_obs_over_n_prior_b = 1 +
      params$i_first_obs_certainty *
        (1 - (i_first_obs_est / pop)),
    hosp_wday_effect_prior_alpha =
      params$hosp_wday_effect_prior_alpha,
    mean_initial_exp_growth_rate_prior_mean =
      params$mean_initial_exp_growth_rate_prior_mean,
    mean_initial_exp_growth_rate_prior_sd =
      params$mean_initial_exp_growth_rate_prior_sd,
    sigma_initial_exp_growth_rate_prior_mode =
      params$sigma_initial_exp_growth_rate_prior_mode,
    sigma_initial_exp_growth_rate_prior_sd =
      params$sigma_initial_exp_growth_rate_prior_sd,
    mode_sigma_ww_site_prior_mode = params$mode_sigma_ww_site_prior_mode,
    mode_sigma_ww_site_prior_sd = params$mode_sigma_ww_site_prior_sd,
    sd_log_sigma_ww_site_prior_mode =
      params$sd_log_sigma_ww_site_prior_mode,
    sd_log_sigma_ww_site_prior_sd =
      params$sd_log_sigma_ww_site_prior_sd,
    eta_sd_sd = params$eta_sd_sd,
    eta_sd_mean = params$eta_sd_mean,
    sigma_i_first_obs_prior_mode = params$sigma_i_first_obs_prior_mode,
    sigma_i_first_obs_prior_sd = params$sigma_i_first_obs_prior_sd,
    p_hosp_prior_mean = params$p_hosp_mean,
    p_hosp_sd_logit = params$p_hosp_sd_logit,
    p_hosp_w_sd_sd = params$p_hosp_w_sd_sd,
    ww_site_mod_sd_sd = params$ww_site_mod_sd_sd,
    inf_feedback_prior_logmean = params$infection_feedback_prior_logmean,
    inf_feedback_prior_logsd = params$infection_feedback_prior_logsd,
    sigma_rt_prior = params$sigma_rt_prior,
    log_phi_g_prior_mean = params$log_phi_g_prior_mean,
    log_phi_g_prior_sd = params$log_phi_g_prior_sd,
    offset_ref_log_r_t_prior_mean = params$offset_ref_log_r_t_prior_mean,
    offset_ref_log_r_t_prior_sd = params$offset_ref_log_r_t_prior_sd,
    offset_ref_logit_i_first_obs_prior_mean =
      params$offset_ref_logit_i_first_obs_prior_mean,
    offset_ref_logit_i_first_obs_prior_sd =
      params$offset_ref_logit_i_first_obs_prior_sd,
    offset_ref_initial_exp_growth_rate_prior_mean =
      params$offset_ref_initial_exp_growth_rate_prior_mean,
    offset_ref_initial_exp_growth_rate_prior_sd =
      params$offset_ref_initial_exp_growth_rate_prior_sd
  )

  return(stan_data_list)
}

#' Get the integer sizes of the wastewater input data
#'
#' @param ww_data Input wastewater dataframe containing one row
#' per observation, with outliers already removed
#' @param lod_col_name  A string representing the name of the
#' column in the input_ww_data that provides a 0 if the data point is not above
#' the LOD and a 1 if the data is below the LOD, default value is `below_LOD`
#'
#' @return A list containing the integer sizes of the follow variables that
#' the stan model requires:
#' owt: number of wastewater observations
#' n_censored: number of censored wastewater observations (below the LOD)
#' n_uncensored: number of uncensored wastewter observations (above the LOD)
#' n_ww_sites: number of wastewater sites
#' n_ww_lab_sites: number of unique wastewater site-lab combinations
#'
#' @export
get_ww_data_sizes <- function(ww_data,
                              lod_col_name = "below_lod") {
  ww_data_present <- nrow(ww_data) != 0
  if (isTRUE(ww_data_present)) {
    # Test for presence of column names
    stopifnot(
      "LOD column name isn't present in input dataset" =
        lod_col_name %in% colnames(ww_data)
    )

    # Number of wastewater observations
    owt <- nrow(ww_data)
    # Number of censored wastewater observations
    n_censored <- sum(ww_data[lod_col_name] == 1)
    # Number of uncensored wastewater observations
    n_uncensored <- owt - n_censored

    # Number of ww sites
    n_ww_sites <- dplyr::n_distinct(ww_data$site_index)

    # Number of unique combinations of wastewater sites and labs
    n_ww_lab_sites <- dplyr::n_distinct(ww_data$lab_site_index)

    data_sizes <- list(
      owt = owt,
      n_censored = n_censored,
      n_uncensored = n_uncensored,
      n_ww_sites = n_ww_sites,
      n_ww_lab_sites = n_ww_lab_sites
    )
  } else {
    data_sizes <- list(
      owt = 0,
      n_censored = 0,
      n_uncensored = 0,
      n_ww_sites = 0,
      n_ww_lab_sites = 0
    )
  }


  return(data_sizes)
}

#' Get wastewater indices and values for stan
#'
#' @param input_ww_data tibble with the input wastewater data and indices
#' needed for stan
#' @param date_time_spine tibble mapping dates to time in days
#' @param lab_site_site_spine tibble mapping lab-sites to sites
#' @param site_subpop_spine tibble mapping sites to subpopulations
#' @param lab_site_subpop_spine tibble mapping lab-sites to subpopulations
#'
#' @return a list of the vectors needed for stan
#' @export
get_ww_indices_and_values <- function(input_ww_data,
                                      date_time_spine,
                                      lab_site_site_spine,
                                      site_subpop_spine,
                                      lab_site_subpop_spine) {
  ww_data_present <- !is.null(input_ww_data)

  # Get a vector of population sizes for each subpop
  subpop_pops <- site_subpop_spine |>
    dplyr::select("subpop_index", "subpop_pop") |>
    dplyr::arrange(.data$subpop_index, "desc") |>
    dplyr::pull(.data$subpop_pop)

  if (isTRUE(ww_data_present)) {
    ww_data_joined <- input_ww_data |>
      dplyr::left_join(date_time_spine, by = "date") |>
      dplyr::left_join(site_subpop_spine, by = c("site_index", "site")) |>
      dplyr::mutate("ind_rel_to_sampled_times" = dplyr::row_number())

    owt <- nrow(ww_data_joined)


    # Get the vector of log LOD values corresponding to each observation
    ww_lod <- ww_data_joined |>
      dplyr::pull("log_lod")


    # Get the vector of log wastewater concentrations
    log_conc <- ww_data_joined |>
      dplyr::pull("log_genome_copies_per_ml")

    # Get censored and uncensored indices, which are relative to the vector
    # of sampled times (e.g. 1:owt)
    ww_censored <- ww_data_joined |>
      dplyr::filter(.data$below_lod == 1) |>
      dplyr::pull(.data$ind_rel_to_sampled_times)
    ww_uncensored <- ww_data_joined |>
      dplyr::filter(.data$below_lod == 0) |>
      dplyr::pull(.data$ind_rel_to_sampled_times)
    stopifnot(
      "Length of censored vectors incorrect" =
        length(ww_censored) + length(ww_uncensored) == owt
    )

    ww_sampled_times <- ww_data_joined |> dplyr::pull("t")
    ww_sampled_subpops <- ww_data_joined |> dplyr::pull("subpop_index")
    lab_site_to_subpop_spine <- lab_site_site_spine |>
      dplyr::left_join(site_subpop_spine, by = "site_index") |>
      pull("subpop_index")
    ww_sampled_lab_sites <- ww_data_joined |> dplyr::pull("lab_site_index")


    ww_values <- list(
      ww_lod = ww_lod,
      subpop_pops = subpop_pops,
      log_conc = log_conc,
      ww_censored = ww_censored,
      ww_uncensored = ww_uncensored,
      ww_sampled_times = ww_sampled_times,
      ww_sampled_subpops = ww_sampled_subpops,
      ww_sampled_lab_sites = ww_sampled_lab_sites
    )
  } else {
    ww_values <- list(
      ww_lod = numeric(),
      subpop_pops = subpop_pops,
      log_conc = numeric(),
      ww_censored = numeric(),
      ww_uncensored = numeric(),
      ww_sampled_times = numeric(),
      ww_sampled_subpops = numeric(),
      ww_sampled_lab_sites = numeric()
    )
  }
  return(ww_values)
}


#' Add time indexing to count data
#'
#' @param input_count_data data frame with dates and counts,
#' but without time indexing.
#'
#' @return The same data frame, with an added
#' time index, including NA rows if dates internal
#' to the timeseries are missing admissions data.
#' @export
#'
#' @examples
#' hosp_data_example <- tibble::tibble(
#'   date = lubridate::ymd("2024-01-01", "2024-01-02", "2024-01-06"),
#'   daily_hosp_admits = c(5, 3, 8)
#' )
#' hosp_data_w_t <- add_time_indexing(hosp_data_example)
add_time_indexing <- function(input_count_data) {
  date_df <- tibble::tibble(date = seq(
    from = min(input_count_data$date),
    to = max(input_count_data$date),
    by = "days"
  )) |>
    dplyr::mutate(t = dplyr::row_number())

  count_data <- input_count_data |>
    dplyr::left_join(date_df, by = "date") |>
    arrange(.data$date)

  return(count_data)
}


#' Get count data integer sizes for stan
#'
#' @param input_count_data a dataframe with the input count data
#' @param forecast_date string indicating the forecast date
#' @param forecast_horizon integer indicating the number of days to make a
#' forecast for
#' @param calibration_time integer indicating the max duration in days that
#' the model is calibrated to hospital admissions for
#' @param last_count_data_date string indicating the date of the last observed
#' count data point
#' @param uot integer indicating the time of model initialization when there are
#' no observations
#' @param count_col_name A string represeting the name of the column in the
#' input_count_data that indicates the number of daily counts,
#' default is `count`
#'
#' @return A list containing the integer sizes of the follow variables that
#' the stan model requires:
#' ht:  integer indicating horizon time for the model(hospital admissions
#' nowcast + forecast time in days)
#' ot: integer indicating the total duration of time that model for producing
#' counts (e.g. cases or admissions) has available calibration data
#' oht: integer indicating the number of count observations
#' n_weeks: number of weeks (rounded up) that counts are generated
#' from the model
#' tot_weeks: number of week(rounded up) that infections are generated for
#' @export
get_count_data_sizes <- function(input_count_data,
                                 forecast_date,
                                 forecast_horizon,
                                 calibration_time,
                                 last_count_data_date,
                                 uot,
                                 count_col_name = "count") {
  nowcast_time <- as.integer(
    lubridate::ymd(forecast_date) - last_count_data_date
  )
  ot <- calibration_time
  ht <- nowcast_time + forecast_horizon
  oht <- input_count_data |>
    dplyr::filter(!is.na(.data[[count_col_name]])) |>
    nrow()
  n_weeks <- ceiling((ot + ht) / 7)
  tot_weeks <- ceiling((ot + uot + ht) / 7)
  count_data_sizes <- list(
    ht = ht,
    ot = calibration_time,
    oht = oht,
    n_weeks = n_weeks,
    tot_weeks = tot_weeks
  )
  return(count_data_sizes)
}
#' Get count data indices
#'
#' @param input_count_data a dataframe with the input count data
#'
#' @return A list containing the vectors of indices that
#' the stan model requires:
#' count_times: a vector of integer times corresponding to the times when the
#' count observations were made
#' @export
get_count_indices <- function(input_count_data) {
  count_times <- input_count_data |>
    dplyr::pull(t)

  count_indices <- list(
    count_times = count_times
  )
  return(count_indices)
}

#' Get count values
#'
#' @param input_count_data a dataframe with the input count data
#' @param ot integer indicating the total duration of time that the
#' model has available calibration data in days
#' @param ht integer indicating the number of days to produce count estimates
#' outside the calibration period (forecast + nowcast time) in days
#' @param count_col_name A string representing the name of the column in the
#' input_count_data that indicates the number of daily counts of the
#' epidemiological indicator, e.g. cases or hospital admissions,
#' default is `count`
#'
#' @return A list containing the necessary vectors of values that
#' the stan model requires:
#' counts: a vector of number of daily count observations
#' day_of_week: a vector indicating the day of the week of each of the dates
#' in the calibration and forecast period
#
#' @export
get_count_values <- function(input_count_data,
                             ot,
                             ht,
                             count_col_name = "count") {
  counts <- input_count_data |>
    dplyr::pull({{ count_col_name }})

  full_dates <- seq(
    from = min(input_count_data$date),
    to = min(input_count_data$date) + lubridate::days(ht + ot - 1),
    by = "days"
  )
  day_of_week <- lubridate::wday(full_dates, week_start = 1)

  count_values <- list(
    counts = counts,
    day_of_week = day_of_week
  )
  return(count_values)
}
