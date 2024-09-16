#' Given a set of prior parameters and stan data, initialize the model
#' near the center of the prior distribution
#'
#' @param stan_data a list of data elements that will be passed to stan
#' @param stdev a numeric value indicating the standard deviation to sample
#' from when initializing, particularly from a standard normal. Also acts as
#' a multiplier on the prior standard deviation, to restrict the initial value
#' to be sampled from the center of the prior. Default is 0.01
#'
#' @return a list of initial values for each of the parameters in the
#' `wwinference` model
get_inits_for_one_chain <- function(stan_data, stdev = 0.01) {
  # Define some variables
  pop <- stan_data$state_pop
  n_weeks <- as.numeric(stan_data$n_weeks)
  tot_weeks <- as.numeric(stan_data$tot_weeks)
  ot <- as.numeric(stan_data$ot)
  ht <- as.numeric(stan_data$ht)
  n_subpops <- as.numeric(stan_data$n_subpops)
  n_ww_lab_sites <- as.numeric(stan_data$n_ww_lab_sites)
  # Estimate of number of initial infections
  i_first_obs_est <- (
    mean(stan_data$hosp[1:7], na.rm = TRUE) / stan_data$p_hosp_prior_mean
  )

  logit_i_frac_est <- stats::qlogis(i_first_obs_est / pop)

  n_subpops <- as.numeric(stan_data$n_subpops)
  n_ww_lab_sites <- as.numeric(stan_data$n_ww_lab_sites)

  init_list <- list(
    w = stats::rnorm(n_weeks - 1, 0, stdev),
    eta_sd = abs(stats::rnorm(1, 0, stdev)),
    eta_i_first_obs = abs(stats::rnorm(n_subpops, 0, stdev)),
    sigma_i_first_obs = abs(stats::rnorm(1, 0, stdev)),
    eta_initial_exp_growth_rate = abs(stats::rnorm(n_subpops, 0, stdev)),
    sigma_initial_exp_growth_rate = abs(stats::rnorm(1, 0, stdev)),
    autoreg_rt = abs(stats::rnorm(
      1,
      stan_data$autoreg_rt_a /
        (stan_data$autoreg_rt_a + stan_data$autoreg_rt_b),
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
    i_first_obs_over_n =
      stats::plogis(stats::rnorm(1, logit_i_frac_est), 0.05),
    mean_initial_exp_growth_rate = stats::rnorm(1, 0, stdev),
    inv_sqrt_phi_h = 1 / sqrt(200) + stats::rnorm(1, 1 / 10000, 1 / 10000),
    mode_sigma_ww_site = abs(stats::rnorm(
      1, stan_data$mode_sigma_ww_site_prior_mode,
      stdev * stan_data$mode_sigma_ww_site_prior_sd
    )),
    sd_log_sigma_ww_site = abs(stats::rnorm(
      1, stan_data$sd_log_sigma_ww_site_prior_mode,
      stdev * stan_data$sd_log_sigma_ww_site_prior_sd
    )),
    eta_log_sigma_ww_site = abs(stats::rnorm(n_ww_lab_sites, 0, stdev)),
    p_hosp_mean = stats::rnorm(
      1, stats::qlogis(stan_data$p_hosp_prior_mean),
      stdev
    ),
    p_hosp_w = stats::rnorm(tot_weeks, 0, stdev),
    p_hosp_w_sd = abs(stats::rnorm(1, 0.01, 0.001)),
    t_peak = stats::rnorm(
      1, stan_data$viral_shedding_pars[1],
      stdev * stan_data$viral_shedding_pars[2]
    ),
    viral_peak = stats::rnorm(
      1, stan_data$viral_shedding_pars[3],
      stdev * stan_data$viral_shedding_pars[4]
    ),
    dur_shed = stats::rnorm(
      1, stan_data$viral_shedding_pars[5],
      stdev * stan_data$viral_shedding_pars[6]
    ),
    log10_g = stats::rnorm(1, stan_data$log10_g_prior_mean, 0.5),
    ww_site_mod_raw = abs(stats::rnorm(n_ww_lab_sites, 0, stdev)),
    ww_site_mod_sd = abs(stats::rnorm(1, 0, stdev)),
    hosp_wday_effect = to_simplex(abs(
      stats::rnorm(7, 1 / 7, stdev)
    )),
    infection_feedback = abs(stats::rnorm(1, 500, 20)),
    # Spatial inits
    log_sigma_generalized = stats::rnorm(1, log(0.01^(n_subpops - 1)), 0.5),
    log_phi = stats::rnorm(1, log(0.05), 0.1),
    log_scaling_factor = stats::rnorm(1, log(1), 0.1),
    non_cent_spatial_dev_ns_mat = matrix(
      stats::rnorm((n_subpops - 1) * n_weeks,
        mean = 0,
        sd = stdev
      ),
      (n_subpops - 1),
      n_weeks
    ),
    norm_vec_aux_site = stats::rnorm(n_weeks, 0, stdev),
    # Initialize the cholesky decomposition matrix if inferring
    # unstructured correlation matrix
    L_Omega = as.matrix(diag(2))
  )

  if (stan_data$corr_structure_switch == 2) {
    init_list$L_Omega <- diag((n_subpops - 1))
  }

  return(init_list)
}
