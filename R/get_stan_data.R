#' Get stan data
#'
#' @param model_type string indicating which model we are getting data for
#' Options are `ww` or `hosp`
#' @param forecast_date string indicating the forecast date
#' @param forecast_time integer indicating the number of days to make a forecast
#' for
#' @param calibration_time integer indicating the max duration in days that
#' the model is calibrated to hospital admissions for
#' @param input_ww_data a dataframe with the input wastewater data
#' @param input_hosp_data a dataframe with the input hospital admissions data
#' @param generation_interval a vector with a zero-truncated normalized pmf of
#' the generation interval
#' @param inf_to_hosp a vector with a normalized pmf of the delay from infection
#'  to hospital admissions
#' @param infection_feedback_pmf a vector with a normalized pmf dictating the
#' delay of infection feedback
#' @param params a dataframe of parameter names and numeric values
#' @param compute_likelihood indicator variable telling stan whether or not to
#' compute the likelihood, default = `1`
#' @param ww_outlier_col_name A string representing the name of the
#' column in the input_ww_data that provides a 0 if the data point is not an
#' outlier to be excluded from the model fit, or a 1 if it is to be excluded
#' default value is `flag_as_ww_outlier`
#' @param lod_col_name A string representing the name of the
#' column in the input_ww_data that provides a 0 if the data point is not above
#' the LOD and a 1 if the data is below the LOD, default value is `below_LOD`
#' @param ww_measurement_col_name A string representing the name of the column
#' in the input_ww_data that indicates the wastewater measurement value in
#' natural scale, default is `ww`
#' @param ww_value_lod_col_name A string representing the name of the column
#' in the input_ww_data that indicates the value of the LOD in natural scale,
#' default is `lod_sewage`
#' @param hosp_value_col_name A string representing the name of the column
#'  in the input_hosp-data that indicates the number of daily hospital
#'  admissions, default is `daily_hosp_admits`
#'
#' @return a list named variables to pass to stan
#' @export
get_stan_data <- function(model_type,
                          forecast_date,
                          forecast_time,
                          calibration_time,
                          input_ww_data,
                          input_hosp_data,
                          generation_interval,
                          inf_to_hosp,
                          infection_feedback_pmf,
                          params,
                          compute_likelihood = 1,
                          ww_outlier_col_name = "flag_as_ww_outlier",
                          lod_col_name = "below_LOD",
                          ww_measurement_col_name = "ww",
                          ww_value_lod_col_name = "lod_sewage",
                          hosp_value_col_name = "daily_hosp_admits") {
  # Assign parameter names
  par_names <- colnames(params)
  for (i in seq_along(par_names)) {
    assign(par_names[i], as.double(params[i]))
  }

  # Indicator variable whether or not to include ww in likelihood
  include_ww <- ifelse(model_type == "ww", 1, 0)

  last_hosp_data_date <- get_last_hosp_data_date(input_hosp_data)

  # Get state pop
  pop <- input_hosp_data |>
    dplyr::select(pop) |>
    unique() |>
    dplyr::pull(pop)

  stopifnot(
    "More than one population size in training data" =
      length(pop) == 1
  )


  if (include_ww == 1) {
    # Test for presence of column names
    stopifnot(
      "Outlier column name isn't present in input dataset" =
        ww_outlier_col_name %in% colnames(input_ww_data)
    )

    # Test to see if ww_data_present
    ww_data_present <- nrow(input_ww_data) != 0
    if (ww_data_present == FALSE) {
      message("No wastewater data present")
    }

    # Filter out wastewater outliers and arrange data for indexing
    ww_data <- input_ww_data |>
      dplyr::filter({{ ww_outlier_col_name }} != 1) |>
      dplyr::arrange(date, site_index)

    ww_data_sizes <- get_ww_data_sizes(
      ww_data,
      lod_col_name
    )
    ww_indices <- get_ww_data_indices(ww_data,
      input_hosp_data,
      owt = ww_data_sizes$owt,
      lod_col_name = lod_col_name
    )
    ww_values <- get_ww_values(
      ww_data,
      ww_measurement_col_name,
      ww_value_lod_col_name,
    )

    stopifnot(
      "Wastewater sampled times not equal to length of input ww data" =
        length(ww_indices$ww_sampled_times) == ww_data_sizes$owt
    )


    message(
      "Prop of population size covered by wastewater: ",
      sum(ww_values$pop_ww) / pop
    )

    # Logic to determine the number of subpopulations to estimate R(t) for:
    # First determine if we need to add an additional subpopulation
    add_auxiliary_site <- ifelse(pop >= sum(ww_values$pop_ww), TRUE, FALSE)
    # Then get the number of subpopulations, the population to normalize by
    # (sum of the subpopulations), and the vector of sizes of each subpopulation
    subpop_data <- get_subpop_data(add_auxiliary_site,
      state_pop = pop,
      pop_ww = ww_values$pop_ww,
      n_ww_sites = ww_data_sizes$n_ww_sites
    )
  } else { # Hospital admissions only model)
    # Still need to specify wastewater input data, so set as 0s. Won't get
    # used by stan to compute the likelihood. None of these will be used.
    owt <- 1
    ww_sampled_times <- c(1)
    log_conc <- c(1)
  }

  # Get the remaining things needed for both models
  hosp_data <- add_time_indexing(input_hosp_data)
  hosp_data_sizes <- get_hosp_data_sizes(
    input_hosp_data = hosp_data,
    forecast_date = forecast_date,
    forecast_time = forecast_time,
    calibration_time = calibration_time,
    last_hosp_data_date = last_hosp_data_date,
    uot = uot,
    hosp_value_col_name = hosp_value_col_name
  )
  hosp_indices <- get_hosp_indices(hosp_data)
  hosp_values <- get_hosp_values(
    hosp_data,
    ot = hosp_data_sizes$ot,
    ht = hosp_data_sizes$ht,
    hosp_value_col_name
  )

  if (include_ww == 1) {
    message(
      "Removed ", nrow(input_ww_data) - ww_data_sizes$owt,
      " outliers from WW data"
    )
  }

  # matrix to transform IHR from weekly to daily
  ind_m <- get_ind_m(
    hosp_data_sizes$ot + hosp_data_sizes$ht,
    hosp_data_sizes$n_weeks
  )
  # matrix to transform p_hosp RW from weekly to daily
  p_hosp_m <- get_ind_m(
    uot + hosp_data_sizes$ot + hosp_data_sizes$ht,
    hosp_data_sizes$tot_weeks
  )

  # Estimate of number of initial infections
  i0 <- mean(hosp_values$hosp_admits[1:7], na.rm = TRUE) / p_hosp_mean

  # package up parameters for stan data object
  viral_shedding_pars <- c(
    t_peak_mean, t_peak_sd, viral_peak_mean, viral_peak_sd,
    duration_shedding_mean, duration_shedding_sd
  )

  hosp_delay_max <- length(inf_to_hosp)

  if (model_type == "ww") {
    data_renewal <- list(
      gt_max = gt_max,
      hosp_delay_max = hosp_delay_max,
      inf_to_hosp = inf_to_hosp,
      dur_inf = dur_inf,
      mwpd = ml_of_ww_per_person_day,
      ot = hosp_data_sizes$ot,
      n_subpops = subpop_data$n_subpops,
      n_ww_sites = ww_data_sizes$n_ww_sites,
      n_ww_lab_sites = ww_data_sizes$n_ww_lab_sites,
      owt = ww_data_sizes$owt,
      oht = hosp_data_sizes$oht,
      n_censored = ww_data_sizes$n_censored,
      n_uncensored = ww_data_sizes$n_uncensored,
      uot = uot,
      ht = hosp_data_sizes$ht,
      n_weeks = hosp_data_sizes$n_weeks,
      ind_m = ind_m,
      tot_weeks = hosp_data_sizes$tot_weeks,
      p_hosp_m = p_hosp_m,
      generation_interval = generation_interval,
      ts = 1:gt_max,
      state_pop = pop,
      subpop_size = subpop_data$subpop_size,
      norm_pop = subpop_data$norm_pop,
      ww_sampled_times = ww_indices$ww_sampled_times,
      hosp_times = hosp_indices$hosp_times,
      ww_sampled_lab_sites = ww_indices$ww_sampled_lab_sites,
      ww_log_lod = ww_values$ww_lod,
      ww_censored = ww_indices$ww_censored,
      ww_uncensored = ww_indices$ww_uncensored,
      hosp = hosp_values$hosp_admits,
      day_of_week = hosp_values$day_of_week,
      log_conc = ww_values$log_conc,
      compute_likelihood = compute_likelihood,
      include_ww = include_ww,
      include_hosp = 1,
      if_l = length(infection_feedback_pmf),
      infection_feedback_pmf = infection_feedback_pmf,
      # All the priors!
      viral_shedding_pars = viral_shedding_pars, # tpeak, viral peak, dur_shed
      autoreg_rt_a = autoreg_rt_a,
      autoreg_rt_b = autoreg_rt_b,
      autoreg_rt_site_a = autoreg_rt_site_a,
      autoreg_rt_site_b = autoreg_rt_site_b,
      autoreg_p_hosp_a = autoreg_p_hosp_a,
      autoreg_p_hosp_b = autoreg_p_hosp_b,
      inv_sqrt_phi_prior_mean = inv_sqrt_phi_prior_mean,
      inv_sqrt_phi_prior_sd = inv_sqrt_phi_prior_sd,
      r_prior_mean = r_prior_mean,
      r_prior_sd = r_prior_sd,
      log10_g_prior_mean = log10_g_prior_mean,
      log10_g_prior_sd = log10_g_prior_sd,
      i0_over_n_prior_a = 1 + i0_certainty * (i0 / pop),
      i0_over_n_prior_b = 1 + i0_certainty * (1 - (i0 / pop)),
      wday_effect_prior_mean = wday_effect_prior_mean,
      wday_effect_prior_sd = wday_effect_prior_sd,
      initial_growth_prior_mean = initial_growth_prior_mean,
      initial_growth_prior_sd = initial_growth_prior_sd,
      sigma_ww_site_prior_mean_mean = sigma_ww_site_prior_mean_mean,
      sigma_ww_site_prior_mean_sd = sigma_ww_site_prior_mean_sd,
      sigma_ww_site_prior_sd_mean = sigma_ww_site_prior_sd_mean,
      sigma_ww_site_prior_sd_sd = sigma_ww_site_prior_sd_sd,
      eta_sd_sd = eta_sd_sd,
      sigma_i0_prior_mode = sigma_i0_prior_mode,
      sigma_i0_prior_sd = sigma_i0_prior_sd,
      p_hosp_prior_mean = p_hosp_mean,
      p_hosp_sd_logit = p_hosp_sd_logit,
      p_hosp_w_sd_sd = p_hosp_w_sd_sd,
      ww_site_mod_sd_sd = ww_site_mod_sd_sd,
      inf_feedback_prior_logmean = infection_feedback_prior_logmean,
      inf_feedback_prior_logsd = infection_feedback_prior_logsd,
      sigma_rt_prior = sigma_rt_prior,
      log_phi_g_prior_mean = log_phi_g_prior_mean,
      log_phi_g_prior_sd = log_phi_g_prior_sd,
      ww_sampled_sites = ww_indices$ww_sampled_sites,
      lab_site_to_site_map = ww_indices$lab_site_to_site_map
    )
  } else if (model_type == "hosp") {
    data_renewal <- list(
      gt_max = gt_max,
      hosp_delay_max = hosp_delay_max,
      inf_to_hosp = inf_to_hosp,
      dur_inf = dur_inf, # this is used bc drift
      mwpd = ml_of_ww_per_person_day,
      ot = hosp_data_sizes$ot,
      owt = owt,
      oht = hosp_data_sizes$oht,
      uot = uot,
      ht = hosp_data_sizes$ht,
      n_weeks = hosp_data_sizes$n_weeks,
      ind_m = ind_m,
      tot_weeks = hosp_data_sizes$tot_weeks,
      p_hosp_m = p_hosp_m,
      generation_interval = generation_interval,
      ts = 1:gt_max,
      n = pop,
      hosp_times = hosp_indices$hosp_times,
      ww_sampled_times = ww_sampled_times,
      hosp = hosp_values$hosp_admits,
      day_of_week = hosp_values$day_of_week,
      log_conc = log_conc,
      compute_likelihood = compute_likelihood,
      include_ww = include_ww,
      include_hosp = 1,
      if_l = length(infection_feedback_pmf),
      infection_feedback_pmf = infection_feedback_pmf,
      # Priors
      viral_shedding_pars = viral_shedding_pars, # tpeak, viral peak,
      # duration shedding
      autoreg_rt_a = autoreg_rt_a,
      autoreg_rt_b = autoreg_rt_b,
      autoreg_p_hosp_a = autoreg_p_hosp_a,
      autoreg_p_hosp_b = autoreg_p_hosp_b,
      inv_sqrt_phi_prior_mean = inv_sqrt_phi_prior_mean,
      inv_sqrt_phi_prior_sd = inv_sqrt_phi_prior_sd,
      r_prior_mean = r_prior_mean,
      r_prior_sd = r_prior_sd,
      log10_g_prior_mean = log10_g_prior_mean,
      log10_g_prior_sd = log10_g_prior_sd,
      i0_over_n_prior_a = 1 + i0_certainty * (i0 / pop),
      i0_over_n_prior_b = 1 + i0_certainty * (1 - (i0 / pop)),
      wday_effect_prior_mean = wday_effect_prior_mean,
      wday_effect_prior_sd = wday_effect_prior_sd,
      initial_growth_prior_mean = initial_growth_prior_mean,
      initial_growth_prior_sd = initial_growth_prior_sd,
      sigma_ww_prior_mean = sigma_ww_site_prior_mean_mean,
      eta_sd_sd = eta_sd_sd,
      p_hosp_prior_mean = p_hosp_mean,
      p_hosp_sd_logit = p_hosp_sd_logit,
      p_hosp_w_sd_sd = p_hosp_w_sd_sd,
      inf_feedback_prior_logmean = infection_feedback_prior_logmean,
      inf_feedback_prior_logsd = infection_feedback_prior_logsd
    )
  } else {
    cli::cli_abort("Unknown model")
    data_renewal <- list()
  }

  stopifnot(
    "Model type not specified properly" =
      !purrr::is_empty(data_renewal)
  )

  return(data_renewal)
}
