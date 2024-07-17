#' Function to generate a global weekly R(t) estimate with added noise
#'
#' @param r_in_weeks The mean R(t) for each week of the simulation, can be
#' longer than the 7*n_weeks
#' @param n_weeks The number of weeks that we want to simulate the data for
#' @param n_days The number of days to generate an R(t) estimate for,
#' can be less than `n_weeks` *7
#' @param global_rt_sd The variation in the R(t) estimate to add
#' intrinsic variability to the infection dynamics
#'
#' @return a weekly R(t) estimate with added noise
get_global_rt <- function(r_in_weeks,
                          n_weeks,
                          global_rt_sd) {
  # Generate the state level weekly R(t) before infection feedback
  # Adds a bit of noise, can add more...
  r_weeks_w_noise <- (r_in_weeks * rnorm(
    length(r_in_weeks),
    1, global_rt_sd
  ))[1:n_weeks]

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
  log_r_site <- matrix(nrow = n_sites + 1, ncol = length(r_weeks))

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
#' @param infection_feedback boolean indicating whether or not to include
#' infection feedback
#' @param infection_feedback_pmf vector of simplex describing the delay from
#' incident infections to feedback on incident infections
#'
#' @return A list containing 3 outputs: i_n: n_subpop X total days matrix of
#' daily incident infections per capita in each subpopulation, r_site: adjusted
#' subpopulation level R(t) estimate, i_n_global: vector of daily incident
#' infections per capita in the global population
#' @export
subpop_inf_process <- function(n_subpops,
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
                               infection_feedback_pmf) {
  i_n <- matrix(nrow = n_subpops, ncol = (ot + ht + uot))
  r_site <- matrix(nrow = n_subpops, ncol = (ot + ht))

  # Generate site level initial infections and growth rates
  initial_growth_site <- rnorm(
    n = n_subpops, mean = initial_growth,
    sd = initial_growth_prior_sd
  )
  log_i0_over_n_site <- rnorm(
    n = n_subpops, mean = log_i0_over_n,
    sd = sd_i0_over_n
  )

  # Set up matrix to convert from weekly to daily
  ind_m <- get_ind_m((ot + ht), n_weeks)
  i_n_global <- rep(0, (uot + ot + ht)) # Global infections
  for (i in 1:(n_subpops)) {
    unadj_r_site_daily <- ind_m %*% unadj_r_site[i, ] # daily R site

    site_output <- model$functions$generate_infections(
      unadj_r_site_daily, # Daily unadjusted R(t) in each site
      uot, # the duration of initialization time for exponential growth
      rev(generation_interval), # the reversed generation interval
      log_i0_over_n_site[i], # log of the initial infections per capita
      initial_growth_site[i], # initial exponential growth rate
      ht, # time after last observed hospital admission
      as.numeric(infection_feedback), # binary indicating if infection feedback
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
#' @export
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
