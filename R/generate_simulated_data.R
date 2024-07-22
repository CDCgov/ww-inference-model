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
#' @param infection_feedback Boolean indicating whether or not to include
#' infection feedback into the infection process, default is `TRUE`
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
                                    mean_log_lod = 3.8,
                                    sd_log_lod = 0.2,
                                    global_rt_sd = 0.03,
                                    sigma_eps = 0.05,
                                    sd_i0_over_n = 0.5,
                                    infection_feedback = TRUE,
                                    input_params_path =
                                      fs::path_package("extdata",
                                        "example_params.toml",
                                        package = "wwinference"
                                      )) {
  # Some quick checks to make sure the inputs work as expected-----------------
  check_rt_length(r_in_weeks, ot, nt, forecast_horizon)
  check_ww_site_pops(pop_size, ww_pop_sites)
  check_site_and_lab_indices(site, lab)

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

  # Get pop fractions of each subpop. There will n_sites + 1 subpops
  pop_fraction <- c(
    ww_pop_sites / pop_size,
    (pop_size - sum(ww_pop_sites)) / pop_size
  )
  n_subpops <- length(pop_fraction)
  # Pull parameter values into memory
  params <- get_params(input_params_path) # load in a single row tibble
  par_names <- colnames(params) # pull them into memory
  for (i in seq_along(par_names)) {
    assign(par_names[i], as.double(params[i]))
  }

  # Create a tibble that maps sites, labs, and population sizes of the sites
  n_sites <- length(unique(site))
  site_lab_map <- data.frame(
    site,
    lab
  ) |>
    dplyr::mutate(
      lab_site = dplyr::row_number()
    ) |>
    dplyr::left_join(data.frame(
      site = 1:n_sites,
      ww_pop = ww_pop_sites
    ))
  n_lab_sites <- nrow(site_lab_map)

  # Define some time variables
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
    dplyr::filter(t == ot + nt) |>
    dplyr::pull(date)

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
    max = gt_max, meanlog = mu_gi, sdlog = sigma_gi, fun_dist = rlnorm, n = 5e6
  ) |> drop_first_and_renormalize()

  # Set infection feedback to generation interval
  infection_feedback_pmf <- generation_interval

  ## Delay from infection to hospital admission ----------------------------
  # incubation period + # time from symptom onset to hospital admission

  # Get incubation period for COVID.
  inc <- make_incubation_period_pmf(
    backward_scale, backward_shape, r
  )
  sym_to_hosp <- make_hospital_onset_delay_pmf(neg_binom_mu, neg_binom_size)

  # Final infection to hospital admissions delay distribution
  inf_to_hosp <- make_reporting_delay_pmf(inc, sym_to_hosp)

  ## Shedding kinetics delay distribution-------------------------------
  vl_trajectory <- model$functions$get_vl_trajectory(
    t_peak_mean, viral_peak_mean,
    duration_shedding_mean, gt_max
  )

  # Global undadjusted R(t) ---------------------------------------------------
  unadj_r_weeks <- get_global_rt(
    r_in_weeks = r_in_weeks,
    n_weeks = n_weeks,
    global_rt_sd = global_rt_sd
  )

  # Subpopulation level R(t)-----------------------------------------------
  # get the matrix of subpop level R(t) estimates, assuming normally distributed
  # around the the log of the state R(t) with stdev of sigma_eps
  unadj_r_site <- subpop_rt_process(
    n_subpops = n_subpops,
    r_weeks = unadj_r_weeks,
    subpop_level_rt_variation = sigma_eps
  )
  # Alternatively, can replace this with
  # r_site <- spatial_rt_process(input_params) #nolint

  # Subpopulation infection dynamics-------------------------------------
  # Function takes in all of the requirements to generation incident infections
  # and R(t) estimates for the unobserved time, calibration, and forecast time
  inf_and_subpop_rt <- subpop_inf_process(
    generate_inf_fxn = model$functions$generate_infections,
    n_subpops = n_subpops,
    uot = uot,
    ot = ot,
    ht = ht,
    unadj_r_site = unadj_r_site,
    initial_growth = initial_growth,
    initial_growth_prior_sd =
      initial_growth_prior_sd,
    i0_over_n = i0_over_n,
    sd_i0_over_n = sd_i0_over_n,
    generation_interval =
      generation_interval,
    infection_feedback =
      infection_feedback,
    infection_feedback_pmf =
      infection_feedback_pmf,
    pop_fraction = pop_fraction
  )

  new_i_over_n_site <- inf_and_subpop_rt$i_n
  r_site <- inf_and_subpop_rt$r_site
  new_i_over_n <- inf_and_subpop_rt$i_n_global



  # Generate expected state level hospitalizations from subpop infections -----

  ## Generate a time varying P(hosp|infection)----------------------------------
  p_hosp_days <- get_time_varying_daily_ihr(
    p_hosp_mean,
    uot,
    ot,
    ht,
    tot_weeks,
    p_hosp_w_sd_sd
  )

  ## Latent per capita admissions--------------------------------------------
  model_hosp_over_n <- model$functions$convolve_dot_product(
    p_hosp_days * new_i_over_n, # individuals who will be hospitalized
    rev(inf_to_hosp),
    uot + ot + ht
  )[(uot + 1):(uot + ot + ht)]


  ## Add weekday effect on hospital admissions-------------------------------
  pred_hosp <- pop_size * model$functions$day_of_week_effect(
    model_hosp_over_n,
    day_of_week_vector,
    hosp_wday_effect
  )
  ## Add observation errror---------------------------------------------------
  # This is negative binomial but could swap out for a different obs error
  pred_obs_hosp <- rnbinom(
    n = length(pred_hosp), mu = pred_hosp,
    size = 1 / ((inv_sqrt_phi_prior_mean)^2)
  )



  # Generate expected observed concentrations from infections in each site-----
  ## Genomes per person per day in each site----------------------------------

  log_g_over_n_site <- get_pred_subpop_gen_per_n(
    convolve_fxn = model$functions$convolve_dot_product,
    n_sites = n_sites,
    uot = uot,
    ot = ot,
    ht = ht,
    new_i_over_n_site = new_i_over_n_site,
    log10_g_prior_mean = log10_g_prior_mean,
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
    ml_of_ww_per_person_day = ml_of_ww_per_person_day
  )

  ## Downsample to simulate reporting/collection process---------------------

  log_obs_conc_lab_site <- downsample_ww_obs(
    log_conc_lab_site = log_conc_lab_site,
    n_lab_sites = n_lab_sites,
    ot = ot,
    ht = ht,
    nt = nt,
    lab_site_reporting_freq = lab_site_reporting_freq,
    lab_site_reporting_latency = lab_site_reporting_latency
  )


  # Global adjusted R(t) --------------------------------------------------
  # I(t)/convolve(I(t), g(t)) #nolint
  # This is not used directly, but we want to have it for comparing to the
  # fit.
  rt <- (new_i_over_n / model$functions$convolve_dot_product(
    new_i_over_n, rev(generation_interval), uot + ot + ht
  ))[(uot + 1):(uot + ot + ht)]


  # Format the data-----------------------------------------------------------

  ww_data <- as.data.frame(t(log_obs_conc_lab_site)) |>
    dplyr::mutate(t = 1:(ot + ht)) |>
    tidyr::pivot_longer(!t,
      names_to = "lab_site",
      names_prefix = "V",
      values_to = "log_conc"
    ) |>
    dplyr::mutate(
      lab_site = as.integer(lab_site)
    ) |>
    dplyr::left_join(date_df, by = "t") |>
    dplyr::left_join(site_lab_map,
      by = "lab_site"
    ) |>
    dplyr::left_join(
      data.frame(
        lab_site = 1:n_lab_sites,
        lod_sewage = lod_lab_site
      ),
      by = c("lab_site")
    ) |> # Remove below LOD values
    dplyr::mutate(
      lod_sewage =
        dplyr::case_when(
          is.na(log_conc) ~ NA,
          !is.na(log_conc) ~ lod_sewage
        )
    ) |>
    dplyr::mutate(
      genome_copies_per_ml = exp(log_conc),
      lod = exp(lod_sewage)
    ) |>
    dplyr::filter(!is.na(genome_copies_per_ml)) |>
    dplyr::rename(site_pop = ww_pop) |>
    dplyr::arrange(site, lab, date) |>
    dplyr::select(date, site, lab, genome_copies_per_ml, lod, site_pop)

  # Make a hospital admissions dataframe for model calibration
  hosp_data <- tibble::tibble(
    t = 1:ot,
    daily_hosp_admits = pred_obs_hosp[1:ot],
    state_pop = pop_size
  ) |>
    dplyr::left_join(
      date_df,
      by = "t"
    ) |>
    dplyr::select(
      date,
      daily_hosp_admits,
      state_pop
    )

  # Make another one for model evaluation
  hosp_data_eval <- tibble::tibble(
    t = 1:(ot + ht),
    daily_hosp_admits_for_eval = pred_obs_hosp,
    state_pop = pop_size
  ) |>
    dplyr::left_join(
      date_df,
      by = "t"
    ) |>
    dplyr::select(
      date,
      daily_hosp_admits_for_eval,
      state_pop
    )

  example_data <- list(
    ww_data = ww_data,
    hosp_data = hosp_data,
    hosp_data_eval = hosp_data_eval
  )

  return(example_data)
}
