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
      date > last_count_data_date - lubridate::days(calibration_time)
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
  ww_data_present <- nrow(preprocessed_ww_data) != 0
  if (ww_data_present == FALSE) {
    message("No wastewater data present")
  }

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
  check_req_ww_cols_present(preprocessed_ww_data,
    conc_col_name = "genome_copies_per_ml",
    lod_col_name = "lod"
  )

  # Filter out wastewater outliers, and remove extra wastewater
  # data. Arrange data for indexing. This is what will be returned.
  ww_data <- preprocessed_ww_data |>
    dplyr::filter(
      exclude != 1,
      date > last_count_data_date -
        lubridate::days(calibration_time)
    ) |>
    dplyr::arrange(date, lab_site_index)

  ww_data_sizes <- get_ww_data_sizes(
    ww_data,
    lod_col_name = "below_lod"
  )

  ww_indices <- get_ww_data_indices(
    ww_data,
    first_count_data_date,
    owt = ww_data_sizes$owt,
    lod_col_name = "below_lod"
  )

  ww_data <- ww_data |>
    dplyr::mutate(
      t = ww_indices$ww_sampled_times
    )

  return(ww_data)
}




#' Get stan data for ww + hosp model
#'
#' @param input_count_data a dataframe of the count data to be passed
#' directly to stan, , must have the following columns: date, count, total_pop
#' @param input_ww_data a dataframe of the wastewater data to be passed
#' directly to stan, must have the following columns: date, site, lab,
#' genome_copies_per_ml, site_pop, below_lod, and exclude
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
#' @param dist_matrix Distance matrix, n_sites x n_sites, passed to a
#' distance-based correlation function for epsilon. If NULL, use an independence
#' correlation function, for current implementation (i.e. all sites' epsilon
#' values are independent and identically distributed) .
#' @param bool_spatial_comp Switch for whether or not infer
#' site-to-site/"spatial" correlation matrix, currently correlation matrix
#' follows exponential correlation structure.
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
#' stan_args <- get_stan_data(
#'   input_count_data_for_stan,
#'   input_ww_data_for_stan,
#'   forecast_date,
#'   forecast_horizon,
#'   calibration_time,
#'   generation_interval,
#'   inf_to_count_delay,
#'   infection_feedback_pmf,
#'   params,
#'   include_ww,
#'   dist_matrix = NULL,
#'   bool_spatial_comp = FALSE
#' )
get_stan_data <- function(input_count_data,
                          input_ww_data,
                          forecast_date,
                          forecast_horizon,
                          calibration_time,
                          generation_interval,
                          inf_to_count_delay,
                          infection_feedback_pmf,
                          params,
                          include_ww,
                          compute_likelihood = 1,
                          dist_matrix,
                          bool_spatial_comp) {
  # Assign parameter names
  par_names <- colnames(params)
  for (i in seq_along(par_names)) {
    assign(par_names[i], as.double(params[i]))
  }

  # Get the last date that there were observations of the epidemiological
  # indicator (aka cases or hospital admissions counts)
  last_count_data_date <- max(input_count_data$date, na.rm = TRUE)
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

  validate_both_datasets(
    input_count_data,
    input_ww_data,
    calibration_time = calibration_time,
    forecast_date = forecast_date
  )


  # Get the total pop, coming from the larger population generating the
  # count data
  pop <- input_count_data |>
    dplyr::select(total_pop) |>
    unique() |>
    dplyr::pull(total_pop)

  assert_single_value(pop,
    arg = "global population",
    add_err_msg = c(
      "More than one global population size",
      "passed to stan"
    )
  )

  last_count_data_date <- max(input_count_data$date, na.rm = TRUE)
  first_count_data_date <- min(input_count_data$date, na.rm = TRUE)
  # Returns a list of the vectors of lod values, the site population sizes in
  # order of the site index, a vector of observations of the log of
  # the genome copies per ml
  ww_values <- get_ww_values(
    input_ww_data,
    ww_measurement_col_name = "genome_copies_per_ml",
    ww_lod_value_col_name = "lod",
    ww_site_pop_col_name = "site_pop"
  )
  # Returns a list with the numbers of elements needed for the stan model
  ww_data_sizes <- get_ww_data_sizes(
    input_ww_data,
    lod_col_name = "below_lod"
  )
  # Returns the vectors of indices you need to map latent variables to
  # observations
  ww_indices <- get_ww_data_indices(
    input_ww_data |> dplyr::select(-t),
    first_count_data_date,
    owt = ww_data_sizes$owt,
    lod_col_name = "below_lod"
  )
  # Ensure that both datasets have overlap with one another, are sufficient
  # in length for the specified calibration time, and have proper time indexing

  stopifnot(
    "Wastewater sampled times not equal to length of input ww data" =
      length(ww_indices$ww_sampled_times) == ww_data_sizes$owt
  )

  message(
    "Prop of population size covered by wastewater: ",
    sum(ww_values$pop_ww) / pop
  )

  if (sum(ww_values$pop_ww) / pop > 1) {
    cli::cli_warn(c(
      "The sum of the wastewater site catchment area populations:",
      "is greater than the global population. While the model supports this",
      "we advise checking your input data to ensure it is specified correctly."
    ))
  }


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

  # Get the sizes of all the elements
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
  i0 <- mean(count_values$count[1:7], na.rm = TRUE) / params$p_hosp_mean

  # package up parameters for stan data object
  viral_shedding_pars <- c(
    params$t_peak_mean, params$t_peak_sd, params$viral_peak_mean,
    params$viral_peak_sd, params$duration_shedding_mean,
    params$duration_shedding_sd
  )

  inf_to_count_delay_max <- length(inf_to_count_delay)

  # If user does / doesn't want spatial comps.
  # We can add an extra step here for when spatial desired and dist_matrix
  #   not given.
  if (bool_spatial_comp) {
    ind_corr_func <- 0L
  } else {
    ind_corr_func <- 1L
    # This dist_matrix will not be used, only needed for stan data specs.
    dist_matrix <- matrix(
      0,
      nrow = subpop_data$n_subpops - 1,
      ncol = subpop_data$n_subpops - 1
    )
  }

  stan_args <- list(
    gt_max = params$gt_max,
    hosp_delay_max = inf_to_count_delay_max,
    inf_to_hosp = inf_to_count_delay,
    mwpd = params$ml_of_ww_per_person_day,
    ot = count_data_sizes$ot,
    n_subpops = subpop_data$n_subpops,
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
    subpop_size = subpop_data$subpop_size,
    norm_pop = subpop_data$norm_pop,
    ww_sampled_times = ww_indices$ww_sampled_times,
    hosp_times = count_indices$count_times,
    ww_sampled_lab_sites = ww_indices$ww_sampled_lab_sites,
    ww_log_lod = ww_values$ww_lod,
    ww_censored = ww_indices$ww_censored,
    ww_uncensored = ww_indices$ww_uncensored,
    hosp = count_values$counts,
    day_of_week = count_values$day_of_week,
    log_conc = ww_values$log_conc,
    compute_likelihood = compute_likelihood,
    include_ww = include_ww,
    include_hosp = 1,
    if_l = length(infection_feedback_pmf),
    infection_feedback_pmf = infection_feedback_pmf,
    # All the priors!
    viral_shedding_pars = viral_shedding_pars, # tpeak, viral peak, dur_shed
    autoreg_rt_a = params$autoreg_rt_a,
    autoreg_rt_b = params$autoreg_rt_b,
    autoreg_rt_site_a = params$autoreg_rt_site_a,
    autoreg_rt_site_b = params$autoreg_rt_site_b,
    autoreg_p_hosp_a = params$autoreg_p_hosp_a,
    autoreg_p_hosp_b = params$autoreg_p_hosp_b,
    inv_sqrt_phi_prior_mean = params$inv_sqrt_phi_prior_mean,
    inv_sqrt_phi_prior_sd = params$inv_sqrt_phi_prior_sd,
    r_prior_mean = params$r_prior_mean,
    r_prior_sd = params$r_prior_sd,
    log10_g_prior_mean = params$log10_g_prior_mean,
    log10_g_prior_sd = params$log10_g_prior_sd,
    i0_over_n_prior_a = 1 + params$i0_certainty * (i0 / pop),
    i0_over_n_prior_b = 1 + params$i0_certainty * (1 - (i0 / pop)),
    hosp_wday_effect_prior_alpha = params$hosp_wday_effect_prior_alpha,
    initial_growth_prior_mean = params$initial_growth_prior_mean,
    initial_growth_prior_sd = params$initial_growth_prior_sd,
    sigma_ww_site_prior_mean_mean = params$sigma_ww_site_prior_mean_mean,
    sigma_ww_site_prior_mean_sd = params$sigma_ww_site_prior_mean_sd,
    sigma_ww_site_prior_sd_mean = params$sigma_ww_site_prior_sd_mean,
    sigma_ww_site_prior_sd_sd = params$sigma_ww_site_prior_sd_sd,
    eta_sd_sd = params$eta_sd_sd,
    sigma_i0_prior_mode = params$sigma_i0_prior_mode,
    sigma_i0_prior_sd = params$sigma_i0_prior_sd,
    p_hosp_prior_mean = params$p_hosp_mean,
    p_hosp_sd_logit = params$p_hosp_sd_logit,
    p_hosp_w_sd_sd = params$p_hosp_w_sd_sd,
    ww_site_mod_sd_sd = params$ww_site_mod_sd_sd,
    inf_feedback_prior_logmean = params$infection_feedback_prior_logmean,
    inf_feedback_prior_logsd = params$infection_feedback_prior_logsd,
    sigma_rt_prior = params$sigma_rt_prior,
    log_phi_g_prior_mean = params$log_phi_g_prior_mean,
    log_phi_g_prior_sd = params$log_phi_g_prior_sd,
    ww_sampled_sites = ww_indices$ww_sampled_sites,
    lab_site_to_site_map = ww_indices$lab_site_to_site_map,
    log_phi_mu_prior = params$log_phi_mu_prior,
    log_phi_sd_prior = params$log_phi_sd_prior,
    l = params$l,
    log_sigma_generalized_mu_prior = params$log_sigma_generalized_mu_prior,
    log_sigma_generalized_sd_prior = params$log_sigma_generalized_sd_prior,
    log_scaling_factor_mu_prior = params$log_scaling_factor_mu_prior,
    log_scaling_factor_sd_prior = params$log_scaling_factor_sd_prior,
    dist_matrix = dist_matrix,
    ind_corr_func = ind_corr_func
  )


  return(stan_args)
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

#' Get wastewater data indices
#'
#' @param ww_data Input wastewater dataframe containing one row
#' per observation, with outliers already removed
#' @param first_count_data_date The earliest day with an observation in the '
#' count dataset, in ISO8601 format YYYY-MM-DD
#' @param owt number of wastewater observations
#' @param lod_col_name A string representing the name of the
#' column in the input_ww_data that provides a 0 if the data point is not above
#' the LOD and a 1 if the data is below the LOD, default value is `below_LOD`
#'
#' @return A list containing the necessary vectors of indices that
#' the stan model requires:
#' ww_censored: the vector of time points that the wastewater observations are
#' censored (below the LOD) in order of the date and the site index
#' ww_uncensored: the vector of time points that the wastewater observations are
#' uncensored (above the LOD) in order of the date and the site index
#' ww_sampled_times: the vector of time points that the wastewater observations
#' are passed in in log_conc in order of the date and the site index
#' ww_sampled_sites: the vector of sites that correspond to the observations
#' passed in in log_conc in order of the date and the site index
#' ww_sampled_lab_sites: the vector of unique combinations of site and labs
#' that correspond to the observations passed in in log_conc in order of the
#' date and the site index
#' lab_site_to_site_map: the vector of sites that correspond to each lab-site
#' @export
get_ww_data_indices <- function(ww_data,
                                first_count_data_date,
                                owt,
                                lod_col_name = "below_lod") {
  # Vector of indices along the list of wastewater concentrations that
  # correspond to censored observations
  ww_data_present <- nrow(ww_data) != 0

  if (isTRUE(ww_data_present)) {
    ww_data_with_index <- ww_data |>
      dplyr::mutate(ind_rel_to_sampled_times = dplyr::row_number())
    ww_censored <- ww_data_with_index |>
      dplyr::filter(.data[[lod_col_name]] == 1) |>
      dplyr::pull(ind_rel_to_sampled_times)
    ww_uncensored <- ww_data_with_index |>
      dplyr::filter(.data[[lod_col_name]] == 0) |>
      dplyr::pull(ind_rel_to_sampled_times)
    stopifnot(
      "Length of censored vectors incorrect" =
        length(ww_censored) + length(ww_uncensored) == owt
    )


    # Need to get the times of wastewater sampling, starting at the first
    # day of hospital admissions data
    ww_date_df <- data.frame(
      date = seq(
        from = first_count_data_date,
        to = max(ww_data$date),
        by = "days"
      ),
      t = 1:(as.integer(max(ww_data$date) - first_count_data_date) + 1)
    )

    # Left join the data mapped to time to the wastewater data
    spine_ww <- ww_data |>
      dplyr::left_join(ww_date_df, by = "date")

    # Pull just the vector of times of wastewater observations
    ww_sampled_times <- spine_ww |>
      dplyr::pull(t)

    # Pull just the indexes of the sites that correspond to the vector of
    # sampled times
    ww_sampled_sites <- ww_data$site_index

    # Pull just the indexes of the lab-sites that correspond to the vector of
    # sampled times
    ww_sampled_lab_sites <- ww_data$lab_site_index

    # Need a vector of indices indicating the site for each lab-site
    lab_site_to_site_map <- ww_data |>
      dplyr::select(lab_site_index, site_index) |>
      dplyr::arrange(lab_site_index, "desc") |>
      dplyr::distinct() |>
      dplyr::pull(site_index)

    ww_data_indices <- list(
      ww_censored = ww_censored,
      ww_uncensored = ww_uncensored,
      ww_sampled_times = ww_sampled_times,
      ww_sampled_sites = ww_sampled_sites,
      ww_sampled_lab_sites = ww_sampled_lab_sites,
      lab_site_to_site_map = lab_site_to_site_map
    )
  } else {
    ww_data_indices <- list(
      ww_censored = c(),
      ww_uncensored = c(),
      ww_sampled_times = c(),
      ww_sampled_sites = c(),
      ww_sampled_lab_sites = c(),
      lab_site_to_site_map = c()
    )
  }


  return(ww_data_indices)
}

#' Get wastewater data values
#'
#' @param ww_data Input wastewater dataframe containing one row
#' per observation, with outliers already removed
#' @param ww_measurement_col_name A string representing the name of the column
#' in the input_ww_data that indicates the wastewater measurement value in
#' natural scale, default is `genome_copies_per_ml`
#' @param ww_lod_value_col_name A string representing the name of the column
#' in the ww_data that indicates the value of the LOD in natural scale,
#' default is `lod`
#' @param ww_site_pop_col_name A string representing the name of the column in
#' the ww_data that indicates the number of people represented by that
#' wastewater catchment, default is `site_pop`
#' @param one_pop_per_site a boolean variable indicating if there should only
#' be on catchment area population per site, default is `TRUE` because this is
#' what the stan model expects
#' @param padding_value an smaller numeric value to add to the the
#' concentration measurements to ensure that log transformation will produce
#' real numbers, default value is `1e-8`
#'
#' @return  A list containing the necessary vectors of values that
#' the stan model requires:
#' ww_lod: a vector of the LODs of the corresponding wastewater measurement
#' pop_ww: a vector of the population sizes of the wastewater catchment areas
#' in order of the sites by site_index
#' log_conc: a vector of the log of the wastewater concentration observation
#' @export
get_ww_values <- function(ww_data,
                          ww_measurement_col_name = "genome_copies_per_ml",
                          ww_lod_value_col_name = "lod",
                          ww_site_pop_col_name = "site_pop",
                          one_pop_per_site = TRUE,
                          padding_value = 1e-8) {
  ww_data_present <- nrow(ww_data) != 0

  if (isTRUE(ww_data_present)) {
    # Get the vector of log LOD values corresponding to each observation
    ww_lod <- ww_data |>
      dplyr::pull({{ ww_lod_value_col_name }}) |>
      log()

    # Get a vector of population sizes
    if (isTRUE(one_pop_per_site)) {
      # Want one population per site during the model calibration period,
      # so just take the average across the populations reported for each
      # observation
      pop_ww <- ww_data |>
        dplyr::select(site_index, {{ ww_site_pop_col_name }}) |>
        dplyr::group_by(site_index) |>
        dplyr::summarise(pop_avg = mean(.data[[ww_site_pop_col_name]])) |>
        dplyr::arrange(site_index, "desc") |>
        dplyr::pull(pop_avg)
    } else {
      # Want a vector of length of the number of observations, corresponding to
      # the population at that time
      pop_ww <- ww_data |>
        dplyr::pull({{ ww_site_pop_col_name }})
    }


    # Get the vector of log wastewater concentrations
    log_conc <- ww_data |>
      dplyr::mutate(
        log_conc =
          (log(.data[[ww_measurement_col_name]] + padding_value))
      ) |>
      dplyr::pull(log_conc)

    ww_values <- list(
      ww_lod = ww_lod,
      pop_ww = pop_ww,
      log_conc = log_conc
    )
  } else {
    ww_values <- list(
      ww_lod = c(),
      pop_ww = c(),
      log_conc = c()
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
    dplyr::arrange(date)

  return(count_data)
}

#' Get subpopulation data
#'
#' @param add_auxiliary_site Boolean indicating whether to add another
#' subpopulation in addition to the wastewater sites to estimate R(t) of
#' @param state_pop The state population size
#' @param pop_ww The population size in each of the wastewater sites
#' @param n_ww_sites The number of wastewater sites
#'
#' @return A list containing the necessary integers and vectors that stan
#' needs to estiamte infection dynamics for each subpopulation
#' @export
#'
#' @examples subpop_data <- get_subpop_data(TRUE, 100000, c(1000, 500), 2)
get_subpop_data <- function(add_auxiliary_site,
                            state_pop,
                            pop_ww,
                            n_ww_sites) {
  if (add_auxiliary_site) {
    # In most cases, wastewater catchment coverage < entire state.
    # So here we add a subpopulation that represents the population not
    # covered by wastewater surveillance
    norm_pop <- state_pop
    n_subpops <- n_ww_sites + 1
    subpop_size <- c(pop_ww, state_pop - sum(pop_ww))
  } else {
    message("Sum of wastewater catchment areas is greater than state pop")
    norm_pop <- sum(pop_ww)
    # If sum catchment areas > state pop,
    # use sum of catchment area pop to normalize
    n_subpops <- n_ww_sites # Only divide the state into n_site subpops
    subpop_size <- pop_ww
  }

  subpop_data <- list(
    norm_pop = norm_pop,
    n_subpops = n_subpops,
    subpop_size = subpop_size
  )
  return(subpop_data)
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
