#' Function to generate a global weekly R(t) estimate with added noise
#'
#' @param r_in_weeks The mean R(t) for each week of the simulation in natural
#' scale, can be longer than the n_weeks
#' @param n_weeks The number of weeks that we want to simulate the data for
#' @param global_rt_sd The variation in the R(t) estimate to add
#' intrinsic variability to the infection dynamics
#' @param zero_padding Numeric value to replace any negative R(t) values with.
#' @return a weekly R(t) estimate with added noise
get_global_rt <- function(r_in_weeks,
                          n_weeks,
                          global_rt_sd,
                          zero_padding = 1e-8) {
  # Generate the state level weekly R(t) before infection feedback
  # Adds a bit of noise, can add more...
  r_weeks_w_noise <- (r_in_weeks * rnorm(
    length(r_in_weeks),
    1, global_rt_sd
  ))[1:n_weeks]
  # Replace negative values if any with very small number
  r_weeks_w_noise[r_weeks_w_noise < 0] <- zero_padding

  return(r_weeks_w_noise)
}


#' Get subpopulation level R(t) estimates assuming time and space independence
#'
#' @param n_subpops integer indicating the number of subpopulations,
#' usually this will be n_sites + 1
#' @param r_weeks The "global" R(t) in weeks
#' @param subpop_level_rt_variation  The standard deviation of the Gaussian
#' to generate deviation at the site level in the R(t) estimate
#'
#' @return A n_sites+1 by n_weeks matrix containing the subpopulation R(t)s
subpop_rt_process <- function(n_subpops,
                              r_weeks,
                              subpop_level_rt_variation) {
  # Initialize output matrix
  log_r_site <- matrix(nrow = n_subpops, ncol = length(r_weeks))

  log_r_global_week <- log(r_weeks)
  for (i in 1:(n_subpops)) {
    # This creates each R(t) vector for each subpopulation, by sampling
    # from a normal distribution centered on the state R(t).
    # In the model, this is an AR(1) process based on the previous deviation
    log_r_site[i, ] <- rnorm(
      n = length(r_weeks),
      mean = log_r_global_week,
      sd = subpop_level_rt_variation
    )
  }
  r_site <- exp(log_r_site)

  return(r_site)
}

#' Get the subpopulation level incident infections
#'
#' @param generate_inf_fxn function indicating how to generate infections,
#' This will typically take the `generate_infections()` function from the
#' `wwinference.stan` model.
#' @param n_subpops integer indicating the number of subpopulations
#' @param uot integer indicating the days for exponential growth initialization
#' to occur (referred to as unobserved time)
#' @param ot integer indicating the number of days we will have observed data
#' for in the calibration period
#' @param ht integer indicating the time after the last observed time to
#' forecast
#' @param unadj_r_site n_subpop x n_weeks matrix of the unadjusted site level
#' R(t)
#' @param initial_growth float indicating the mean of the initial growth rate
#'  in the unobserved time
#' @param initial_growth_prior_sd float indicating the standard deviation on
#' the initial growth rate across subpopulations
#' @param i0_over_n float indicating the mean of the initial per capita
#' infections
#' @param sd_i0_over_n float indicating the standard deviation of the log of the
#' initial infections per capita across subpopulations
#' @param generation_interval vector of simplex describing the probability of
#' each time from infection to onwards transmission
#' @param infection_feedback numeric indicating the strength of the infection
#' feedback
#' @param infection_feedback_pmf vector of simplex describing the delay from
#' incident infections to feedback on incident infections
#' @param pop_fraction vector of a simplex of length n_subpops, indicating
#' the proportion of the global population that subpopulation represents
#'
#' @return A list containing 3 outputs: i_n: n_subpop X total days matrix of
#' daily incident infections per capita in each subpopulation, r_site: adjusted
#' subpopulation level R(t) estimate, i_n_global: vector of daily incident
#' infections per capita in the global population
#'
subpop_inf_process <- function(generate_inf_fxn,
                               n_subpops,
                               uot,
                               ot,
                               ht,
                               unadj_r_site,
                               initial_growth,
                               initial_growth_prior_sd,
                               i0_over_n,
                               sd_i0_over_n,
                               generation_interval,
                               infection_feedback,
                               infection_feedback_pmf,
                               pop_fraction) {
  i_n <- matrix(nrow = n_subpops, ncol = (ot + ht + uot))
  r_site <- matrix(nrow = n_subpops, ncol = (ot + ht))

  # Generate site level initial infections and growth rates
  initial_growth_site <- rnorm(
    n = n_subpops, mean = initial_growth,
    sd = initial_growth_prior_sd
  )
  # This is initial infections at the first observation time
  log_i0_over_n_site <- rnorm(
    n = n_subpops, mean = log(i0_over_n),
    sd = sd_i0_over_n
  )


  n_weeks <- ncol(unadj_r_site)
  # Set up matrix to convert from weekly to daily
  ind_m <- get_ind_m((ot + ht), n_weeks)
  i_n_global <- rep(0, (uot + ot + ht)) # Global infections
  for (i in 1:(n_subpops)) {
    unadj_r_site_daily <- ind_m %*% unadj_r_site[i, ] # daily R site

    site_output <- generate_inf_fxn(
      unadj_r_site_daily, # Daily unadjusted R(t) in each site
      uot, # the duration of initialization time for exponential growth
      rev(generation_interval), # the reversed generation interval
      log_i0_over_n_site[i], # log of the initial infections per capita
      initial_growth_site[i], # initial exponential growth rate
      ht, # time after last observed hospital admission
      as.numeric(infection_feedback), # strength of infection feedback
      rev(infection_feedback_pmf) # reversed inf feedback delay pmf
    )
    # matrix to hold infections
    i_n[i, ] <- site_output[[1]]

    # Adjusted R(t) estimate in each site
    r_site[i, ] <- site_output[[2]]

    # Cumulatively sum infections to get overall state infections
    i_n_global <- i_n_global + pop_fraction[i] * site_output[[1]]
  }

  output <- list(
    i_n = i_n,
    r_site = r_site,
    i_n_global = i_n_global
  )
  return(output)
}

#' Get time varying IHR
#'
#' @param p_hosp_mean mean probability of hospital admission given infection
#' @param uot integer indicating the days for exponential growth initialization
#' to occur (referred to as unobserved time)
#' @param ot integer indicating the number of days we will have observed data
#' for in the calibration period
#' @param ht integer indicating the time after the last observed time to
#' forecast
#' @param tot_weeks integer indicating the total number of weeks including the
#' `uot`
#' @param p_hosp_w_sd_sd the standard deviation in logit scale of the
#' IHR
#'
#' @return a vector of daily IHR values
get_time_varying_daily_ihr <- function(p_hosp_mean,
                                       uot,
                                       ot,
                                       ht,
                                       tot_weeks,
                                       p_hosp_w_sd_sd) {
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

  return(p_hosp_days)
}

#' Get the predicted genomes per person in each subpopulation
#'
#' @param convolve_fxn function used to convolve infections with delay pmf.
#' This will typically take `convolve_dot_product()` function from the
#' `wwinference.stan` model
#' @param n_sites integer indicating the number of sites
#' @param uot integer indicating the days for exponential growth initialization
#' to occur (referred to as unobserved time)
#' @param ot integer indicating the number of days we will have observed data
#' for in the calibration period
#' @param ht integer indicating the time after the last observed time to
#' @param new_i_over_n_site matrix of n_subpops rows x n total time points
#' columns indicating the per capita infections in each subpopulation
#' @param log10_g_prior_mean log10 scale prior on the number of genomes shed
#' per infection
#' @param vl_trajectory vector of simplex indicating the proportion of all
#' shedding occuring on each day after infection
#'
#' @return A matrix of `n_sites` rows and `ot`+ `ht` columns indicating
#' the predicted number of genomes per person per day in each site each day
get_pred_subpop_gen_per_n <- function(convolve_fxn,
                                      n_sites,
                                      uot,
                                      ot,
                                      ht,
                                      new_i_over_n_site,
                                      log10_g_prior_mean,
                                      vl_trajectory) {
  log_g_over_n_site <- matrix(nrow = n_sites, ncol = (ot + ht))

  for (i in 1:n_sites) {
    # Convolve infections with shedding kinetics
    model_net_i <- convolve_fxn(
      new_i_over_n_site[i, ],
      rev(vl_trajectory),
      (uot + ot + ht)
    )[(uot + 1):(uot + ot + ht)]
    # Scale by average genomes shed per infection
    log_g_over_n_site[i, ] <- log(10) * log10_g_prior_mean +
      log(model_net_i + 1e-8)
  }

  return(log_g_over_n_site)
}

#' Get the predicted concentrations in each lab site
#'
#' @param n_lab_sites Integer indicating the number of unique combinations of
#' labs and sites in the dataset.
#' @param ot integer indicating the number of days we will have observed data
#' for in the calibration period
#' @param ht integer indicating the time after the last observed time to
#' @param log_g_over_n_site matrix of n_site rows and ot + ht columns indicating
#' the genomes per person each day in each site
#' @param log_m_lab_sites vector of the lab-site mutlipliers
#' @param sigma_ww_lab_site vector of the lab_site observation errors
#' @param site vector of integers indicating which site (WWTP) each separate
#' lab-site observation comes from
#' @param ml_of_ww_per_person_day Scalar indicating the number of mL of
#' wastewater produced per person per day
#'
#' @return A matrix of `n_lab_sites` rows and `ot` + `ht` columns indcating the
#' predicted concentration of wastewater in each lab-site and day
get_pred_obs_conc <- function(n_lab_sites,
                              ot,
                              ht,
                              log_g_over_n_site,
                              log_m_lab_sites,
                              sigma_ww_lab_site,
                              site,
                              ml_of_ww_per_person_day) {
  log_conc_lab_site <- matrix(nrow = n_lab_sites, ncol = (ot + ht))
  for (i in 1:n_lab_sites) {
    log_g_w_multiplier <- log_g_over_n_site[site[i], ] +
      log_m_lab_sites[i] # Add site level multiplier in log scale
    log_conc_lab_site[i, ] <- log_g_w_multiplier +
      rnorm(
        n = (ot + ht), mean = 0,
        sd = sigma_ww_lab_site[i]
      ) - log(ml_of_ww_per_person_day)
  }

  return(log_conc_lab_site)
}


#' Downsample the predicted wastewater concentrations based on the
#' lab site reporting frequency
#'
#' @param log_conc_lab_site The matrix of n_lab_sites by n time points
#' indicating the underlying expected  observed concentrations
#' @param n_lab_sites Integer indicating the number of unique lab-site
#' combinations
#' @param ot integer indicating the number of days we will have observed data
#' for in the calibration period
#' @param ht integer indicating the time after the last observed time to
#' the end of the forecast time
#' @param nt integer indicating the time after the last observed epi indicator
#'  and before the forecast date, of which there can still be wastewater
#'  observations
#' @param lab_site_reporting_freq vector indicating the mean frequency of
#' wastewater measurements in each site per day (e.g. 1/7 is once per week)

#' @return A sparse matrix of `n_lab_sites` rows and `ot` + `ht` columns of
#' but with NAs for when observations are not measured/reported.
downsample_for_frequency <- function(log_conc_lab_site,
                                     n_lab_sites,
                                     ot,
                                     ht,
                                     nt,
                                     lab_site_reporting_freq) {
  log_obs_conc_lab_site <- matrix(nrow = n_lab_sites, ncol = ot + ht)
  for (i in 1:n_lab_sites) {
    # Get the indices where we observe concentrations
    st <- sample(1:(ot + nt), round((ot + nt) * lab_site_reporting_freq[i]))
    # Calculate log concentration for the days that we have observations
    log_obs_conc_lab_site[i, st] <- log_conc_lab_site[i, st]
  }

  return(log_obs_conc_lab_site)
}

#' Truncate the predicted wastewater concentrations based on the
#' lab site reporting latency and the observed time and horizon time
#'
#' @param log_conc_lab_site The matrix of n_lab_sites by n time points
#' indicating the underlying expected  observed concentrations
#' @param n_lab_sites Integer indicating the number of unique lab-site
#' combinations
#' @param ot integer indicating the number of days we will have observed data
#' for in the calibration period
#' @param ht integer indicating the time after the last observed time to
#' the end of the forecast time
#' @param nt integer indicating the time after the last observed epi indicator
#'  and before the forecast date, of which there can still be wastewater
#'  observations
#' @param lab_site_reporting_latency vector indicating the number of days
#' from the forecast date of the last possible observation

#' @return A sparse matrix of `n_lab_sites` rows and `ot` + `ht` columns of
#' but with NAs for when observations are not measured/reported.
truncate_for_latency <- function(log_conc_lab_site,
                                 n_lab_sites,
                                 ot,
                                 ht,
                                 nt,
                                 lab_site_reporting_latency) {
  log_obs_conc_lab_site <- log_conc_lab_site
  for (i in 1:n_lab_sites) {
    # Get the last day there can be none NAs
    last_index_day <- ot + nt - lab_site_reporting_latency[i]
    # Replace with NAs behond last index day
    log_obs_conc_lab_site[i, last_index_day:(ot + ht)] <- NA
  }

  return(log_obs_conc_lab_site)
}


#' Format the wastewater data as a tidy data frame
#'
#' @param log_obs_conc_lab_site matrix of numeric values where rows are the
#' site-lab combinations and columns are the observed time points
#' @param ot integer indicating the number of days we will have observed data
#' for in the calibration period
#' @param ht integer indicating the time after the last observed epi indicator
#'  and before the forecast date, of which there can still be wastewater
#'  observations
#' @param date_df tibble of columns `date` and `t` that map time in days to
#' dates
#' @param site_lab_map tibble mapping sites, labs, lab site indices, and the
#' population size of the site
#' @param lod_lab_site vector of numerics indicating the LOD in each lab and
#' site combination
#'
#' @return a tidy dataframe containing observed wastewater concentrations
#' in log estimated genome copies per mL for each site and lab at each time
#' point
format_ww_data <- function(log_obs_conc_lab_site,
                           ot,
                           ht,
                           date_df,
                           site_lab_map,
                           lod_lab_site) {
  n_lab_sites <- nrow(site_lab_map)
  ww_data <- as.data.frame(t(log_obs_conc_lab_site)) |>
    dplyr::mutate(t = 1:(ot + ht)) |>
    tidyr::pivot_longer(!t,
      names_to = "lab_site",
      names_prefix = "V",
      values_to = "log_conc"
    ) |>
    dplyr::mutate(
      lab_site = as.integer(.data$lab_site)
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
    ) |>
    dplyr::rename(
      "log_lod" = "lod_sewage",
      "log_genome_copies_per_ml" = "log_conc",
      "site_pop" = "ww_pop"
    ) |>
    # Remove missing values
    dplyr::filter(!is.na(.data$log_genome_copies_per_ml)) |>
    dplyr::arrange(.data$site, .data$lab, .data$date) |>
    dplyr::select(
      "date", "site", "lab", "log_genome_copies_per_ml",
      "log_lod", "site_pop"
    )

  return(ww_data)
}


#' Format the hospital admissions data into a tidy dataframe
#'
#' @param pred_obs_hosp vector of non-negative integers indicating the number
#' of hospital admissions on each day
#' @param dur_obs integer indicating the number of days we want the
#' observations for
#' @param pop_size integer indicating the population size of the admissions
#' catchment area
#' @param date_df tibble of columns `date` and `t` that map time in days to
#' dates
#'
#' @return a tidy dataframe containing counts of admissions by date alongside
#' population size
format_hosp_data <- function(pred_obs_hosp,
                             dur_obs,
                             pop_size,
                             date_df) {
  hosp_data <- tibble::tibble(
    t = 1:dur_obs,
    daily_hosp_admits = pred_obs_hosp[1:dur_obs],
    state_pop = pop_size
  ) |>
    dplyr::left_join(
      date_df,
      by = "t"
    ) |>
    dplyr::select(
      "date",
      "daily_hosp_admits",
      "state_pop"
    )
  return(hosp_data)
}


#' Format the subpopulation-level hospital admissions data into a tidy
#' dataframe
#'
#' @param pred_obs_hosp_subpop matrix of non-negative integers indicating the
#' number of hospital admissions on each day in each subpopulation. Rows are
#' subpopulations, columns are time points
#' @param dur_obs integer indicating the number of days we want the
#' observations for
#' @param subpop_map tibble mapping the numbered subpopulations to the
#' wastewater sites, must contain columns "subpop_index" and "subpop_name"
#' @param date_df tibble of columns `date` and `t` that map time in days to
#' dates
#'
#' @return a tidy dataframe containing counts of admissions by date alongside
#' population size for each subpopulation
format_subpop_hosp_data <- function(pred_obs_hosp_subpop,
                                    dur_obs,
                                    subpop_map,
                                    date_df) {
  subpop_hosp_data <- as.data.frame(t(pred_obs_hosp_subpop)) |>
    dplyr::mutate(t = seq_len(ncol(pred_obs_hosp_subpop))) |>
    dplyr::filter(t <= dur_obs) |>
    tidyr::pivot_longer(!t,
      names_to = "subpop_index",
      names_prefix = "V",
      values_to = "daily_hosp_admits"
    ) |>
    dplyr::left_join(
      date_df,
      by = "t"
    ) |>
    dplyr::left_join(
      subpop_map,
      by = "subpop_index"
    ) |>
    dplyr::select(
      "date",
      "subpop_name",
      "daily_hosp_admits",
      "subpop_pop"
    )
  return(subpop_hosp_data)
}


#' Back- calculate R(t) from incident infections and the generation interval
#'
#' @description
#' Note that the forward renewal process is not a simple convolution of
#' incident infections and the generation interval -- because it assumes that
#' the generation interval being passed in is indexed starting at day 1, and
#' that on any particular index day there is no contribution from individuals
#' infected on that day. This is a reasonable assumption, but to align the
#' implementation of the forward process with our backward process, we have to
#' add a 0 density to day 0 of the passed in generation interval.
#'
#'
#' @param new_i vector of numerics that spans the length of `tot_time`,
#' representing the new incident infections per day
#' @param convolve_fxn function used to convolve infections with delay pmf
#' This will typically take the `convolve_dot_product()` function from the
#' `wwinference.stan` model
#' @param generation_interval vector of simplex describing the probability of
#' each time from infection to onwards transmission
#' @param uot integer indicating the days for exponential growth initialization
#' to occur (referred to as unobserved time)
#' @param tot_time integer indicating the total time we have incident
#' infections for
#'
#' @return a numeric vector of length(`tot_time` - `uot`) that represents
#' the effective reproductive number
#' @export
calc_rt <- function(new_i,
                    convolve_fxn,
                    generation_interval,
                    uot,
                    tot_time) {
  rt <- (new_i / convolve_fxn(
    new_i, rev(c(0, generation_interval)), tot_time # Assert no contributions
    # from those infected on the day of interest we're calculating.
  ))[(uot + 1):tot_time]
  return(rt)
}

#' Create a mapping of sites, labs, and population sizes in each site
#'
#' @param site vector of integers or characters uniquely identifying a site,
#'  or the population in the wastewater catchment area
#' @param lab vector of integer or character uniquely identifying a lab where
#' a sample was processed
#' @param ww_pop_sites vector of integers indicating the population size in
#' each site
#'
#' @return a tibble mapping sites, labs and population sizes
create_site_lab_map <- function(site,
                                lab,
                                ww_pop_sites) {
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

  return(site_lab_map)
}
