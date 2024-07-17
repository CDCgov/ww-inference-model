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
#' @param n_sites integer indicating the number of wastewater sites, assumes
#' we will always have n_sites + 1 subpopulations
#' @param r_weeks The "global" R(t) in weeks
#' @param subpop_level_rt_variation  The standard deviation of the Gaussian
#' to generate deviation at the site level in the R(t) estimate
#'
#' @return A n_sites+1 by n_weeks matrix containing the subpopulation R(t)s
subpop_rt_process <- function(n_sites,
                              r_weeks,
                              subpop_level_rt_variation) {
  # Initialize output matrix
  log_r_site <- matrix(nrow = n_sites + 1, ncol = length(r_weeks))

  log_r_global_week <- log(r_weeks)
  for (i in 1:(n_sites + 1)) {
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
