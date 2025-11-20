#' Generate simulated data from the underlying model's generative process
#' @description
#' Function that allows the user to generate hospital admissions and site-level
#' wastewater data directly from the generative model, specifying the conditions
#' and parameters to generate from.
#'
#' @param r_in_weeks vector indcating the mean weekly R(t) that drives infection
#'  dynamics at the state-level. This gets jittered with random noise to add
#'  week-to-week variation.
#' @param n_sites integer indicating the number of sites
#' @param ww_pop_sites vector indicating the population size in the
#' catchment area in each of those sites (order must match)
#' @param pop_size integer indicating the population size in the hypothetical
#' state, default is `3e6`
#' @param site vector of integers indicating which site (WWTP) each separate
#' lab-site observation comes from
#' @param lab vector of integers indicating which lab the lab-site observations
#' come from
#' @param ot integer indicating the observed time: length of hospital admissions
#'  calibration time in days
#' @param nt integer indicating the nowcast time: length of time between last
#' hospital admissions date and forecast date in days
#' @param forecast_horizon integer indicating the duration of the forecast in
#' days e.g. 28 days
#' @param sim_start_date character string in ISO8601 format YYYY-MM-DD
#'  indicating the start date of the simulation, used to get a weekday vector
#' @param hosp_wday_effect a vector that is a simplex of length 7 describing
#' how the hospital admissions are spread out over a week, starting at
#' Monday = 1
#' @param i0_over_n float between 0 and 1 indicating the initial per capita
#' infections in the state
#' @param initial_growth float indicating the exponential growth rate in
#' infections (daily) during the unobserved time
#' @param sd_in_lab_level_multiplier float indicating the standard deviation in
#' the log of the site-lab level multiplier determining how much variation
#' there is systematically in site-labs from the state mean
#' @param mean_obs_error_in_ww_lab_site float indicating the mean day-to-day
#' variation in observed wastewater concentrations across all lab-sites
#' @param mean_reporting_freq float indicating the mean frequency of wastewater
#'  measurements across sites in per day (e.g. 1/7 is once per week)
#' @param sd_reporting_freq float indicating the standard deviation in the
#' frequency of wastewater measurements across sites
#' @param mean_reporting_latency float indicating the mean time from forecast
#' date to last wastewater sample collection date, across sites
#' @param sd_reporting_latency float indicating the standard deviation in the
#' time from the forecast date to the last wastewater sample collection date,
#' across sites
#' @param mean_log_lod float indicating the mean log of the LOD in each lab-site
#' @param sd_log_lod float indicating the standard deviation in the log of the
#' LOD across sites
#' @param global_rt_sd float indicating the ammount of standard deviation to
#' add to the passed in weekly R(t) to add variability. Default is `0.03`
#' @param sigma_eps float indicating the standard deviation between the log of
#' the state R(t) and the log of the subpopulation R(t) across time, in log
#' scale. Default is `0.05`
#' @param sd_i0_over_n float indicating the standard deviation between log of
#' initial infections per capita, default is `0.5`
#' @param if_feedback Boolean indicating whether or not to include
#' infection feedback into the infection process, default is `FALSE`, which
#' sets the strength of the infection feedback to 0.
#' If `TRUE`, this will apply an infection feedback drawn from the prior.
#' @param subpop_phi Vector of numeric values indicating the overdispersion
#' parameter phi in the hospital admissions observation process in each
#' subpopulation
#' @param input_params_path path to the toml file with the parameters to use
#' to generate the simulated data
#'
#' @return a list containing three dataframes. hosp_data is a dataframe
#' containing the number of daily hospital admissions by day for a theoretical
#' US state, for the duration of the specified calibration period.
#' hosp_data_eval is a dataframe containing the number of daily hospital
#' admissions by day for a theoretical US state, for the entire evaluation
#' period.
#' ww_data is a dataframe containing the measured wastewater concentrations
#' in each site alongside other metadata necessary for modeling that data.
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate a simulated dataset from a hypothetical state with 6 sites and 2
#' # different labs
#' sim_data <- generate_simulated_data(
#'   n_sites = 6,
#'   site = c(1, 2, 3, 4, 5, 6, 6),
#'   lab = c(1, 1, 1, 1, 2, 2, 3),
#'   ww_pop_sites = c(1e5, 4e5, 2e5, 1.5e5, 5e4, 3e5),
#'   pop_size = 2e6
#' )
#' hosp_data <- sim_data$hosp_data
#' ww_data <- sim_data$ww_data
#' }
generate_simulated_data <- function(r_in_weeks = # nolint
                                      c(
                                        rep(1.1, 5), rep(0.9, 5),
                                        1 + 0.007 * 1:16
                                      ),
                                    n_sites = 4,
                                    ww_pop_sites = c(4e5, 2e5, 1e5, 5e4),
                                    pop_size = 3e6,
                                    site = c(1, 1, 2, 3, 4),
                                    lab = c(1, 2, 3, 3, 3),
                                    ot = 90,
                                    nt = 9,
                                    forecast_horizon = 28,
                                    sim_start_date = lubridate::ymd(
                                      "2023-09-01"
                                    ),
                                    hosp_wday_effect = c(
                                      0.95, 1.01, 1.02,
                                      1.02, 1.01, 1,
                                      0.99
                                    ) / 7,
                                    i0_over_n = 5e-4,
                                    initial_growth = 1e-4,
                                    sd_in_lab_level_multiplier = 0.25,
                                    mean_obs_error_in_ww_lab_site = 0.2,
                                    mean_reporting_freq = 1 / 5,
                                    sd_reporting_freq = 1 / 20,
                                    mean_reporting_latency = 7,
                                    sd_reporting_latency = 3,
                                    mean_log_lod = 5,
                                    sd_log_lod = 0.2,
                                    global_rt_sd = 0.03,
                                    sigma_eps = 0.05,
                                    sd_i0_over_n = 0.5,
                                    if_feedback = FALSE,
                                    subpop_phi = c(25, 50, 70, 40, 100),
                                    input_params_path =
                                      fs::path_package("extdata",
                                        "example_params.toml",
                                        package = "wwinference"
                                      )) {
  # Some quick checks to make sure the inputs work as expected-----------------
  assert_rt_correct_length(r_in_weeks, ot, nt, forecast_horizon)
  assert_ww_site_pops_lt_total(pop_size, ww_pop_sites)
  assert_site_lab_indices_align(site, lab)

  # Expose the stan functions into the global environment--------------------
  model <- cmdstanr::cmdstan_model(
    stan_file = system.file(
      "stan", "wwinference.stan",
      package = "wwinference"
    ),
    compile = TRUE,
    compile_standalone = TRUE,
    force_recompile = TRUE
  )
  model$expose_functions(global = TRUE)

  # Get other variables needed for forward simulation ------------------------
  params <- get_params(input_params_path) # load in parameters

  # Get pop fractions of each subpop. There will n_sites + 1 subpops
  pop_fraction <- c(
    ww_pop_sites / pop_size,
    (pop_size - sum(ww_pop_sites)) / pop_size
  )
  n_subpops <- length(pop_fraction)

  # Create a tibble that maps sites, labs, and population sizes of the sites
  n_sites <- length(unique(site))
  site_lab_map <- create_site_lab_map(
    site = site,
    lab = lab,
    ww_pop_sites = ww_pop_sites
  )


  n_lab_sites <- nrow(site_lab_map)

  # Define some time variables
  uot <- params$uot # Time period of exponential growth outside of
  # admissions calibratin window (initialization of renewal process)
  ht <- nt + forecast_horizon
  n_weeks <- ceiling((ot + ht) / 7) # calibration + forecast time
  tot_weeks <- ceiling((uot + ot + ht) / 7) # initialization time +
  # calibration + forecast time


  # We need dates to get a weekday vector
  dates <- seq(
    from = sim_start_date, to =
      (sim_start_date + lubridate::days(ot + nt + ht - 1)), by = "days"
  )
  log_i0_over_n <- log(i0_over_n)
  day_of_week_vector <- lubridate::wday(dates, week_start = 1)
  date_df <- data.frame(
    t = 1:(ot + nt + ht),
    date = dates
  )

  forecast_date <- date_df |>
    dplyr::filter(.data$t == !!ot + !!nt) |>
    dplyr::pull("date")

  # Set the lab-site multiplier presumably from lab measurement processes
  log_m_lab_sites <- rnorm(n_lab_sites,
    mean = 0, sd = sd_in_lab_level_multiplier
  ) # This is the magnitude shift (multiplier in natural scale) on the
  # observations, presumably from things like concentration method, PCR type,
  # collection type, etc.

  # Assign a site level observation error to each site, but have it scale
  # inversely with the catchment area of the site for now. Eventually, we will
  # want to impose the expected variability as a function of the contributing
  # infections, but since this module isn't currently in the model we will
  # just do this for now.
  sigma_ww_lab_site <- mean(site_lab_map$ww_pop) *
    mean_obs_error_in_ww_lab_site / site_lab_map$ww_pop

  # Set randomly the lab-site reporting avg frequency (per day) and the
  # reporting latency (in days). Will use this to sample times in the observed
  # data
  lab_site_reporting_freq <- abs(rnorm(
    n = n_lab_sites, mean = mean_reporting_freq,
    sd = sd_reporting_freq
  ))
  lab_site_reporting_latency <- pmax(1, ceiling(rnorm(
    n = n_lab_sites,
    mean = mean_reporting_latency, sd = sd_reporting_latency
  )))
  # Set a lab-site-specific LOD in log scale
  lod_lab_site <- rnorm(n_lab_sites, mean = mean_log_lod, sd = sd_log_lod)

  # Delay distributions ------------------------------------------------------
  # Note, these are COVID specific, and we will
  # generally want to specify these outside model configuration

  ## Generation interval------------------------------------------------------
  # Double censored and zero-truncated
  generation_interval <- simulate_double_censored_pmf(
    max = params$gt_max, meanlog = params$mu_gi, sdlog = params$sigma_gi,
    fun_dist = rlnorm, n = 5e6
  ) |> drop_first_and_renormalize()

  # Set infection feedback to generation interval
  infection_feedback_pmf <- generation_interval

  ## Delay from infection to hospital admission ----------------------------
  # incubation period + # time from symptom onset to hospital admission

  # Get incubation period for COVID.
  inc <- make_incubation_period_pmf(
    params$backward_scale, params$backward_shape, params$r
  )
  sym_to_hosp <- make_hospital_onset_delay_pmf(
    params$neg_binom_mu,
    params$neg_binom_size
  )

  # Final infection to hospital admissions delay distribution
  inf_to_hosp <- make_reporting_delay_pmf(inc, sym_to_hosp)

  ## Shedding kinetics delay distribution-------------------------------
  vl_trajectory <- model$functions$get_vl_trajectory(
    params$t_peak_mean, params$viral_peak_mean,
    params$duration_shedding_mean, params$gt_max
  )

  # Global undadjusted R(t) ---------------------------------------------------
  unadj_r_weeks <- get_global_rt(
    r_in_weeks = r_in_weeks,
    n_weeks = n_weeks,
    global_rt_sd = global_rt_sd
  )

  # Daily unadjusted R(t)
  ind_m <- get_ind_m((ot + ht), n_weeks)
  unadj_r_daily <- ind_m %*% unadj_r_weeks

  # Subpopulation level R(t)-----------------------------------------------
  # get the matrix of subpop level R(t) estimates, assuming normally distributed
  # around the the log of the state R(t) with stdev of sigma_eps
  unadj_r_site <- subpop_rt_process(
    n_subpops = n_subpops,
    r_weeks = unadj_r_weeks,
    subpop_level_rt_variation = sigma_eps
  )


  # Subpopulation infection dynamics-------------------------------------
  # Function takes in all of the requirements to generation incident infections
  # and R(t) estimates for the unobserved time, calibration, and forecast time
  if (isTRUE(if_feedback)) {
    infection_feedback <- rlnorm(1,
      meanlog = params$infection_feedback_prior_logmean,
      sdlog = params$infection_feedback_prior_logsd
    )
  } else {
    infection_feedback <- 0
  }

  inf_and_subpop_rt <- subpop_inf_process(
    generate_inf_fxn = model$functions$generate_infections,
    n_subpops = n_subpops,
    uot = uot,
    ot = ot,
    ht = ht,
    unadj_r_site = unadj_r_site,
    initial_growth = initial_growth,
    initial_growth_prior_sd = params$mean_initial_exp_growth_rate_prior_sd,
    i0_over_n = i0_over_n,
    sd_i0_over_n = sd_i0_over_n,
    generation_interval = generation_interval,
    infection_feedback = infection_feedback,
    infection_feedback_pmf = infection_feedback_pmf,
    pop_fraction = pop_fraction
  )

  new_i_over_n_site <- inf_and_subpop_rt$i_n
  r_site <- inf_and_subpop_rt$r_site
  new_i_over_n <- inf_and_subpop_rt$i_n_global


  # Generate expected state level hospitalizations from subpop infections -----

  ## Generate a time varying P(hosp|infection)----------------------------------
  p_hosp_days <- get_time_varying_daily_ihr(
    p_hosp_mean = params$p_hosp_mean,
    uot = params$uot,
    ot = ot,
    ht = ht,
    tot_weeks = tot_weeks,
    p_hosp_w_sd_sd = params$p_hosp_w_sd_sd
  )

  ## Latent per capita admissions--------------------------------------------
  # This won't be used other than for the unit test
  model_hosp_over_n <- model$functions$convolve_dot_product(
    p_hosp_days * new_i_over_n, # individuals who will be hospitalized
    rev(inf_to_hosp),
    uot + ot + ht
  )[(uot + 1):(uot + ot + ht)]

  # Also compute per capita hosps for each subpopulation
  model_hosp_subpop_over_n <- matrix(
    nrow = n_subpops,
    ncol = (ot + ht)
  )
  for (i in 1:n_subpops) {
    model_hosp_subpop_over_n[i, ] <- model$functions$convolve_dot_product(
      p_hosp_days * new_i_over_n_site[i, ],
      rev(inf_to_hosp),
      uot + ot + ht
    )[(uot + 1):(uot + ot + ht)]
  }

  # unit test to make sure these are equivalent
  if (!all.equal(
    colSums(model_hosp_subpop_over_n * pop_fraction),
    model_hosp_over_n,
    tolerance = 1e-8
  )) {
    cli::cli_abort("Sum of convolutions not equal to convolution of sums")
  }


  ## Add weekday effect on hospital admissions-------------------------------
  pred_hosp <- pop_size * model$functions$day_of_week_effect(
    model_hosp_over_n,
    day_of_week_vector,
    hosp_wday_effect
  )

  pred_hosp_subpop <- matrix(
    nrow = n_subpops,
    ncol = (ot + ht)
  )
  for (i in 1:n_subpops) {
    pred_hosp_subpop[i, ] <- pop_fraction[i] * pop_size *
      model$functions$day_of_week_effect(
        model_hosp_subpop_over_n[i, ],
        day_of_week_vector,
        hosp_wday_effect
      )
  }


  ## Add observation error---------------------------------------------------
  # Use negative binomial but could swap out for a different obs error.
  # Each subpopulation has its own dispersion parameter, then we sum
  # the observations to get the population total
  pred_obs_hosp_subpop <- matrix(
    nrow = n_subpops,
    ncol = (ot + ht)
  )
  for (i in 1:n_subpops) {
    pred_obs_hosp_subpop[i, ] <- rnbinom(
      n = length(pred_hosp_subpop[i, ]), mu = pred_hosp_subpop[i, ],
      size = subpop_phi[i]
    )
  }
  pred_obs_hosp <- colSums(pred_obs_hosp_subpop)


  # Generate expected observed concentrations from infections in each site-----
  ## Genomes per person per day in each site----------------------------------

  log_g_over_n_site <- get_pred_subpop_gen_per_n(
    convolve_fxn = model$functions$convolve_dot_product,
    n_sites = n_sites,
    uot = uot,
    ot = ot,
    ht = ht,
    new_i_over_n_site = new_i_over_n_site,
    log10_g_prior_mean = params$log10_g_prior_mean,
    vl_trajectory = vl_trajectory
  )
  ## Site-lab level observation error----------------------------------------
  log_conc_lab_site <- get_pred_obs_conc(
    n_lab_sites = n_lab_sites,
    ot = ot,
    ht = ht,
    log_g_over_n_site = log_g_over_n_site,
    log_m_lab_sites = log_m_lab_sites,
    sigma_ww_lab_site = sigma_ww_lab_site,
    site = site,
    ml_of_ww_per_person_day = params$ml_of_ww_per_person_day
  )

  ## Downsample to simulate reporting/collection process---------------------

  # Create evaluation data with same reporting freq but go through the entire
  # time period
  log_obs_conc_lab_site_eval <- downsample_for_frequency(
    log_conc_lab_site = log_conc_lab_site,
    n_lab_sites = n_lab_sites,
    ht = ht,
    ot = ot,
    nt = nt,
    lab_site_reporting_freq = lab_site_reporting_freq
  )


  log_obs_conc_lab_site <- truncate_for_latency(
    log_conc_lab_site = log_obs_conc_lab_site_eval,
    n_lab_sites = n_lab_sites,
    ot = ot,
    ht = ht,
    nt = nt,
    lab_site_reporting_latency = lab_site_reporting_latency
  )


  # Global adjusted R(t) --------------------------------------------------
  # I(t)/convolve(I(t), g(t)) #nolint
  # This is not used directly, but we want to have it for comparing to the
  # fit.
  rt <- calc_rt(
    new_i = new_i_over_n,
    convolve_fxn = model$functions$convolve_dot_product,
    generation_interval = generation_interval,
    uot = uot,
    tot_time = (uot + ot + ht)
  )

  # Format the data-----------------------------------------------------------

  ww_data <- format_ww_data(
    log_obs_conc_lab_site = log_obs_conc_lab_site,
    ot = ot,
    ht = ht,
    date_df = date_df,
    site_lab_map = site_lab_map,
    lod_lab_site = lod_lab_site
  )

  ww_data_eval <- format_ww_data(
    log_obs_conc_lab_site = log_obs_conc_lab_site_eval,
    ot = ot + ht,
    ht = 0,
    date_df = date_df,
    site_lab_map = site_lab_map,
    lod_lab_site = lod_lab_site
  ) |>
    dplyr::rename(
      "log_genome_copies_per_ml_eval" = "log_genome_copies_per_ml"
    )

  # Artificially add values below the LOD----------------------------------
  # Replace it with an NA, will be used as an example of how to format data
  # properly.
  min_ww_val <- min(ww_data$log_genome_copies_per_ml)
  ww_data <- ww_data |>
    dplyr::mutate(
      "log_genome_copies_per_ml" =
        dplyr::case_when(
          .data$log_genome_copies_per_ml ==
            !!min_ww_val ~ 0.5 * .data$log_lod,
          TRUE ~ .data$log_genome_copies_per_ml
        )
    )
  ww_data_eval <- ww_data_eval |>
    dplyr::mutate(
      "log_genome_copies_per_ml_eval" =
        dplyr::case_when(
          .data$log_genome_copies_per_ml_eval ==
            !!min_ww_val ~ 0.5 * .data$log_lod,
          TRUE ~ .data$log_genome_copies_per_ml_eval
        )
    )


  # Make a hospital admissions dataframe for model calibration
  hosp_data <- format_hosp_data(
    pred_obs_hosp = pred_obs_hosp,
    dur_obs = ot,
    pop_size = pop_size,
    date_df = date_df
  )

  hosp_data_eval <- format_hosp_data(
    pred_obs_hosp = pred_obs_hosp,
    dur_obs = (ot + ht),
    pop_size = pop_size,
    date_df = date_df
  ) |>
    dplyr::rename(
      "daily_hosp_admits_for_eval" = "daily_hosp_admits"
    )

  # Make a subpopulation level hospital admissions data
  # For now this will only be used for evaluation, eventually, can add
  # feature to use this in calibration
  subpop_map <- tibble::tibble(
    subpop_index = as.character(1:n_subpops),
    subpop_pop = pop_size * pop_fraction,
    subpop_name = c(1:n_sites, NA)
  ) |>
    dplyr::mutate(subpop_name = ifelse(!is.na(subpop_name),
      glue::glue("Site: {subpop_name}"),
      "remainder of population"
    ))

  subpop_hosp_data <- format_subpop_hosp_data(
    pred_obs_hosp_subpop = pred_obs_hosp_subpop,
    dur_obs = ot,
    subpop_map = subpop_map,
    date_df = date_df
  )

  subpop_hosp_data_eval <- format_subpop_hosp_data(
    pred_obs_hosp_subpop = pred_obs_hosp_subpop,
    dur_obs = (ot + ht),
    subpop_map = subpop_map,
    date_df = date_df
  ) |>
    dplyr::rename(
      "daily_hosp_admits_for_eval" = "daily_hosp_admits"
    )

  # Global R(t)
  true_rt <- tibble::tibble(
    unadj_rt_daily = as.numeric(unadj_r_daily),
    realized_rt = rt
  ) |>
    dplyr::mutate(
      t = 1:(ot + ht)
    ) |>
    dplyr::left_join(date_df, by = "t")


  example_data <- list(
    ww_data = ww_data,
    ww_data_eval = ww_data_eval,
    hosp_data = hosp_data,
    hosp_data_eval = hosp_data_eval,
    subpop_hosp_data = subpop_hosp_data,
    subpop_hosp_data_eval = subpop_hosp_data_eval,
    true_global_rt = true_rt
  )

  return(example_data)
}
