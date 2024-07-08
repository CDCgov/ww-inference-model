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
#' state
#' @param site vector of integers indicating which site (WWTP) each separate
#' lab-site observation comes frm
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
#'   site = c(rep(1, 4), rep(2, 2))
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
                                    pop_size = 1e6,
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
                                    input_params_path =
                                      fs::path_package("extdata",
                                        "example_params.toml",
                                        package = "wwinference"
                                      )) {
  # Some quick checks to make sure the inputs work as expected
  stopifnot(
    "weekly R(t) passed in isn't long enough" =
      length(r_in_weeks) >= (ot + nt + forecast_horizon) / 7
  )
  stopifnot(
    "Sum of wastewater site populations is greater than state pop" =
      pop_size > sum(ww_pop_sites)
  )
  stopifnot(
    "Site and lab indices don't align" =
      length(site) == length(lab)
  )


  # Get pop fractions of each subpop. There will n_sites + 1 subpops
  pop_fraction <- c(
    ww_pop_sites / pop_size,
    (pop_size - sum(ww_pop_sites)) / pop_size
  )

  # Expose the stan functions into the global environment
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

  # Delay distributions: Note, these are COVID specific, and we will
  # generally want to specify these outside model configuration

  # Double censored, zero truncated, generation interval
  generation_interval <- simulate_double_censored_pmf(
    max = gt_max, meanlog = mu_gi, sdlog = sigma_gi, fun_dist = rlnorm, n = 5e6
  ) |> drop_first_and_renormalize()

  # Set infection feedback to generation interval
  infection_feedback_pmf <- generation_interval
  infection_feedback_rev_pmf <- rev(infection_feedback_pmf)
  infection_feedback <- 0
  if_feedback <- 1

  # Delay from infection to hospital admission: incubation period +
  # time from symptom onset to hospital admission

  # Get incubation period for COVID.
  inc <- make_incubation_period_pmf(
    backward_scale, backward_shape, r
  )
  sym_to_hosp <- make_hospital_onset_delay_pmf(neg_binom_mu, neg_binom_size)

  # Infection to hospital admissions delay distribution
  inf_to_hosp <- make_reporting_delay_pmf(inc, sym_to_hosp)

  # Shedding kinetics delay distribution
  vl_trajectory <- model$functions$get_vl_trajectory(
    t_peak_mean, viral_peak_mean,
    duration_shedding_mean, gt_max
  )

  # Generate the state level weekly R(t) before infection feedback
  # Adds a bit of noise, can add more...
  unadj_r_weeks <- (r_in_weeks * rnorm(length(r_in_weeks), 1, 0.03))[1:n_weeks]

  # Convert to daily for input into renewal equation
  ind_m <- get_ind_m(ot + ht, n_weeks)
  unadj_r <- ind_m %*% unadj_r_weeks


  # Generate the site-level expected observed concentrations -----------------
  # first by adding variation to the site-level R(t) in each site,
  # and then adding lab-site level multiplier and observation error


  ### Generate the site level infection dynamics-------------------------------
  new_i_over_n_site <- matrix(nrow = n_sites + 1, ncol = (uot + ot + ht))
  r_site <- matrix(nrow = n_sites + 1, ncol = (ot + ht))
  # Generate site-level R(t)
  log_r_state_week <- log(unadj_r_weeks)
  log_r_site <- matrix(nrow = n_sites + 1, ncol = n_weeks)
  initial_growth_site <- vector(length = n_sites + 1)
  log_i0_over_n_site <- vector(length = n_sites + 1)
  for (i in 1:(n_sites + 1)) {
    # This creates each R(t) vector for each subpopulation, by sampling
    # from a normal distribution centered on the state R(t).
    # In the model, this is an AR(1) process based on the previous deviation
    log_r_site[i, ] <- rnorm(
      n = n_weeks,
      mean = log_r_state_week,
      sd = 0.05
    ) # sigma_rt

    # Generate deviations in the initial growth rate and initial incidence
    initial_growth_site[i] <- rnorm(
      n = 1, mean = initial_growth,
      sd = initial_growth_prior_sd
    )
    # This is I0/N at the first unobserved time
    log_i0_over_n_site[i] <- rnorm(
      n = 1, mean = log_i0_over_n,
      sd = 0.5
    )
  }


  new_i_over_n <- rep(0, (uot + ot + ht)) # State infections
  for (i in 1:(n_sites + 1)) {
    unadj_r_site <- ind_m %*% exp(log_r_site[i, ]) # daily R site

    site_output <- model$functions$generate_infections(
      unadj_r_site, # Daily unadjusted R(t) in each site
      uot, # the duration of initialization time for exponential growth
      rev(generation_interval), # the reversed generation interval
      log_i0_over_n_site[i], # log of the initial infections per capita
      initial_growth_site[i], # initial exponential growth rate
      ht, # time after last observed hospital admission
      infection_feedback, # binary indicating whether or not inf feedback
      infection_feedback_rev_pmf # inf feedback delay pmf
    )
    # matrix to hold infections
    new_i_over_n_site[i, ] <- site_output[[1]]
    # Cumulatively sum infections to get overall state infections
    new_i_over_n <- new_i_over_n + pop_fraction[i] * site_output[[1]]
    # Adjusted R(t) estimate in each site
    r_site[i, ] <- site_output[[2]]
  }

  # State adjusted R(t) = I(t)/convolve(I(t), g(t))
  rt <- (new_i_over_n / model$functions$convolve_dot_product(
    new_i_over_n, rev(generation_interval), uot + ot + ht
  ))[(uot + 1):(uot + ot + ht)]


  # Generate expected state level hospitalizations from subpop infections -----

  # Generate a time varying P(hosp|infection),
  p_hosp_int_logit <- qlogis(p_hosp_mean) # p_hosp_mean is in linear scale
  p_hosp_m <- get_ind_m(uot + ot + ht, tot_weeks) # matrix needed to convert
  # from weekly to daily
  p_hosp_w_logit <- p_hosp_int_logit + rnorm(
    tot_weeks - 1, 0,
    p_hosp_w_sd_sd
  )
  # random walk on p_hosp in logit scale
  p_hosp_logit_weeks <- c(
    p_hosp_int_logit,
    p_hosp_w_logit
  ) # combine with intercept
  p_hosp_logit_days <- p_hosp_m %*% c(
    p_hosp_int_logit,
    p_hosp_w_logit
  ) # convert to days
  p_hosp_days <- plogis(p_hosp_logit_days) # convert back to linear scale


  # Get expected trajectory of hospital admissions from incident infections
  # by convolving scaled incident infections with delay from infection to
  # hospital admission
  model_hosp_over_n <- model$functions$convolve_dot_product(
    p_hosp_days * new_i_over_n, # individuals who will be hospitalized
    rev(inf_to_hosp),
    uot + ot + ht
  )[(uot + 1):(uot + ot + ht)]
  # only care about hospital admission in observed time, but need uot infections

  exp_hosp <- pop_size * model$functions$day_of_week_effect(
    model_hosp_over_n,
    day_of_week_vector,
    hosp_wday_effect
  )
  # Add observation error, get hospital admissions in the forecast period
  exp_obs_hosp <- rnbinom(
    n = length(exp_hosp), mu = exp_hosp,
    size = 1 / ((inv_sqrt_phi_prior_mean)^2)
  )



  ## Generate site-level mean genomes from infections in each site-------
  log_g_over_n_site <- matrix(nrow = n_sites, ncol = (ot + ht))

  for (i in 1:n_sites) {
    # Convolve infections with shedding kinetics
    model_net_i <- model$functions$convolve_dot_product(
      new_i_over_n_site[i, ],
      rev(vl_trajectory),
      (uot + ot + ht)
    )[(uot + 1):(uot + ot + ht)]
    # Scale by average genomes shed per infection
    log_g_over_n_site[i, ] <- log(10) * log10_g_prior_mean +
      log(model_net_i + 1e-8)
  }


  # Add on site-lab-level observation error -----------------------------------
  log_obs_g_over_n_lab_site <- matrix(nrow = n_lab_sites, ncol = (ot + ht))
  for (i in 1:n_lab_sites) {
    log_g_w_multiplier <- log_g_over_n_site[site[i], ] +
      log_m_lab_sites[i] # Add site level multiplier in log scale
    log_obs_g_over_n_lab_site[i, ] <- log_g_w_multiplier +
      rnorm(
        n = (ot + ht), mean = 0,
        sd = sigma_ww_lab_site[i]
      ) # + add observation error in log scale
  }

  # Sample from some lab-sites more frequently than others and add different
  # latencies for each lab-site
  log_obs_conc_lab_site <- matrix(nrow = n_lab_sites, ncol = ot + ht)
  for (i in 1:n_lab_sites) {
    # Get the indices where we observe concentrations
    st <- sample(1:(ot + nt), round((ot + nt) * lab_site_reporting_freq[i]))
    # cut off end based on latency
    stl <- pmin((ot + nt - lab_site_reporting_latency[i]), st)
    # Calculate log concentration for the days that we have observations
    log_obs_conc_lab_site[i, stl] <- log_obs_g_over_n_lab_site[i, stl] -
      log(ml_of_ww_per_person_day)
  }

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
    daily_hosp_admits = exp_obs_hosp[1:ot],
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
    daily_hosp_admits_for_eval = exp_obs_hosp,
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
