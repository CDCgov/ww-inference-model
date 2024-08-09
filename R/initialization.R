#' Given a set of prior parameters and stan data, initialize the model
#' near the center of the prior distribution
#'
#' @param stan_data a list of data elements that will be passed to stan
#' @param params a dataframe of parameter values that are passed to stan
#' to specify the priors in the model
#' @param stdev a numeric value indicating the standard deviation to sample
#' from when initializing, particularly from a standard normal. Also acts as
#' a multiplier on the prior standard deviation, to restrict the initial value
#' to be sampled from the center of the prior. Default is 0.01
#'
#' @return a list of initial values for each of the parameters in the
#' `wwinference` model
get_inits <- function(stan_data, params, stdev = 0.01) {
  # Define some variables
  pop <- stan_data$state_pop
  n_weeks <- as.numeric(stan_data$n_weeks)
  tot_weeks <- as.numeric(stan_data$tot_weeks)
  ot <- as.numeric(stan_data$ot)
  ht <- as.numeric(stan_data$ht)
  n_subpops <- as.numeric(stan_data$n_subpops)
  n_ww_lab_sites <- as.numeric(stan_data$n_ww_lab_sites)
  # Estimate of number of initial infections
  i0 <- mean(stan_data$hosp[1:7], na.rm = TRUE) / params$p_hosp_mean

  n_subpops <- as.numeric(stan_data$n_subpops)
  n_ww_lab_sites <- as.numeric(stan_data$n_ww_lab_sites)

  init_list <- list(
    w = stats::rnorm(n_weeks - 1, 0, stdev),
    eta_sd = abs(stats::rnorm(1, 0, stdev)),
    eta_i0 = abs(stats::rnorm(n_subpops, 0, stdev)),
    sigma_i0 = abs(stats::rnorm(1, 0, stdev)),
    eta_growth = abs(stats::rnorm(n_subpops, 0, stdev)),
    sigma_growth = abs(stats::rnorm(1, 0, stdev)),
    autoreg_rt = abs(stats::rnorm(
      1,
      params$autoreg_rt_a / (params$autoreg_rt_a + params$autoreg_rt_b),
      0.05
    )),
    log_r_mu_intercept = stats::rnorm(
      1,
      convert_to_logmean(1, stdev),
      convert_to_logsd(1, stdev)
    ),
    error_site = matrix(
      stats::rnorm(n_subpops * n_weeks,
        mean = 0,
        sd = stdev
      ),
      n_subpops,
      n_weeks
    ),
    autoreg_rt_site = abs(stats::rnorm(1, 0.5, 0.05)),
    autoreg_p_hosp = abs(stats::rnorm(1, 1 / 100, 0.001)),
    sigma_rt = abs(stats::rnorm(1, 0, stdev)),
    i0_over_n = stats::plogis(stats::rnorm(1, stats::qlogis(i0 / pop), 0.05)),
    initial_growth = stats::rnorm(1, 0, stdev),
    inv_sqrt_phi_h = 1 / sqrt(200) + stats::rnorm(1, 1 / 10000, 1 / 10000),
    sigma_ww_site_mean = abs(stats::rnorm(
      1, params$sigma_ww_site_prior_mean_mean,
      stdev * params$sigma_ww_site_prior_mean_sd
    )),
    sigma_ww_site_sd = abs(stats::rnorm(
      1, params$sigma_ww_site_prior_sd_mean,
      stdev * params$sigma_ww_site_prior_sd_sd
    )),
    sigma_ww_site_raw = abs(stats::rnorm(n_ww_lab_sites, 0, stdev)),
    p_hosp_mean = stats::rnorm(1, stats::qlogis(params$p_hosp_mean), stdev),
    p_hosp_w = stats::rnorm(tot_weeks, 0, stdev),
    p_hosp_w_sd = abs(stats::rnorm(1, 0.01, 0.001)),
    t_peak = stats::rnorm(1, params$t_peak_mean, stdev * params$t_peak_sd),
    viral_peak = stats::rnorm(
      1, params$viral_peak_mean,
      stdev * params$viral_peak_sd
    ),
    dur_shed = stats::rnorm(
      1, params$duration_shedding_mean,
      stdev * params$duration_shedding_sd
    ),
    log10_g = stats::rnorm(1, params$log10_g_prior_mean, 0.5),
    ww_site_mod_raw = abs(stats::rnorm(n_ww_lab_sites, 0, stdev)),
    ww_site_mod_sd = abs(stats::rnorm(1, 0, stdev)),
    hosp_wday_effect = to_simplex(stats::rnorm(7, 1 / 7, stdev)),
    infection_feedback = abs(stats::rnorm(1, 500, 20))
  )
  return(init_list)
}
