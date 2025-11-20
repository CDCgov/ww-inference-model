#' Simulate daily double censored PMF. From epinowcast:
#' https://package.epinowcast.org/dev/reference/simulate_double_censored_pmf.html #nolint
#'
#' This function simulates the probability mass function of a  daily
#' double-censored process. The process involves two distributions: a primary
#' distribution which represents the censoring process for the primary event
#' and another distribution (which is offset by the primary).
#'
#' Based off of:
#' https://www.medrxiv.org/content/10.1101/2024.01.12.24301247v1
#'
#' @param max Maximum value for the computed CDF. If not specified, the maximum
#' value is the maximum simulated delay.
#' @param fun_primary Primary distribution function (default is \code{runif}).
#' @param fun_dist Distribution function to be added to the primary (default is
#' \code{rlnorm}).
#' @param n Number of simulations (default is 1e6).
#' @param primary_args List of additional arguments to be passed to the primary
#' distribution function.
#' @param dist_args List of additional arguments to be passed to the
#' distribution function.
#' @param ... Additional arguments to be passed to the distribution function.
#' This is an alternative to `dist_args`.
#'
#' @return A numeric vector representing the PMF.
simulate_double_censored_pmf <- function(
  max, fun_primary = stats::runif, primary_args = list(),
  fun_dist = stats::rlnorm,
  dist_args = list(...), n = 1e6, ...
) {
  primary <- do.call(fun_primary, c(list(n), primary_args))
  secondary <- primary + do.call(fun_dist, c(list(n), dist_args))
  delay <- floor(secondary) - floor(primary)
  if (missing(max)) {
    max <- base::max(delay)
  }
  cdf <- ecdf(delay)(0:max)
  pmf <- c(cdf[1], diff(cdf))
  vec_outside_tol <- abs(sum(pmf) - 1L) > 1e-10
  while (vec_outside_tol) {
    pmf <- pmf / sum(pmf)
    vec_outside_tol <- abs(sum(pmf) - 1L) > 1e-10
  }
  return(pmf)
}

#' Drop the first element of a simplex and re-normalize the result to sum to 1.
#'
#' When this vector corresponds to the generation interval distribution, we
#' want to drop this first bin. The renewal equation assumes that same-day
#' infection and onward transmission does not occur, and we assume
#' everything is 1 indexed not 0 indeced. We need to
#' manually drop the first element from the PMF vector.
#'
#' @param x A numeric vector, sums to 1. Corresponds to a discretized PDF or PMF
#'   (usually the GI distribution).
#'
#' @return A numeric vector, sums to 1.
drop_first_and_renormalize <- function(x) {
  # Check input sums to 1
  stopifnot(abs(sum(x) - 1) < 1e-8)
  # Drop and renormalize
  y <- x[2:length(x)] / sum(x[2:length(x)])
  vec_outside_tol <- abs(sum(y) - 1L) > 1e-10
  # Normalize until within tolerance
  while (vec_outside_tol) {
    y <- y / sum(y)
    vec_outside_tol <- abs(sum(y) - 1L) > 1e-10
  }
  return(y)
}

#' @title Make incubation period pmf
#' @description When the default arguments are used, this returns a pmf
#' corresponding to the incubation period for COVID after Omicron used in
#' Park et al 2023. These estimates are from early Omicron.
#' @param backward_scale numeric indicating the scale parameter for the Weibull
#' used in producing the incubation period distribution. default is `3.60` for
#' COVID
#' @param backward_shape numeric indicating the shape parameter for the Weibull
#' used in producing the incubation period distribution, default is `1.50` for
#' COVID
#' @param r numeric indicating the exponential rate used in producing the
#' correction on the incubaion period distribution, default is `0.15` for COVID
#'
#' @return pmf of incubation period
make_incubation_period_pmf <- function(backward_scale = 3.60,
                                       backward_shape = 1.50,
                                       r = 0.15) {
  # From: Park, Sang Woo, et al. "Inferring the differences in incubation-period
  # and generation-interval distributions of the Delta and Omicron variants of
  # SARS-CoV-2." Proceedings of the National Academy of Sciences 120.22 (2023):
  # e2221887120.

  # "However, when we account for growth-rate differences and reestimate the
  # forward incubation periods, we find that both variants have similar
  # incubation-period distributions with a mean of 4.1 d (95% CI: 3.8 to 4.4 d)
  # for the Delta variant and 4.2 d (95% CI: 3.6 to 4.9 d) for the Omicron
  # variant Fig. 3B)."

  # Fits a Weibull to the data

  # Relies on fundamental assumption about epidemic growth rate.


  discr_gr_adj_weibull <- tibble::tibble(
    time = seq(0, 23, by = 1), # 23 seems to get most of the distribution mass
    density0 = dweibull(time,
      shape = backward_shape,
      scale = backward_scale
    ) * exp(r * time)
  )

  inc_period_pmf <- wwinference::to_simplex(discr_gr_adj_weibull$density0)
  return(inc_period_pmf)
}


#' @title Make hospital onset delay pmf
#' @description Uses the parameter estimates from cfa-parameter-estimates,
#' which is based on Danache et al linelist data from symptom onset to hospital
#' admission. See below:
#' https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0261428
#'
#' @param neg_binom_mu float indicating the mean of the negative binomial shaped
#' delay from symptom onset to hospital admissions, default is `6.98665` from
#' fit to data in above paper
#' @param neg_binom_size float indicating the dispersion parameter in the
#' negative binomial delay from symptom onset to hospital admissions, default
#' is `2.490848` from fit to data in above paper
#'
#' @return pmf of distribution from symptom onset to hospital admission
make_hospital_onset_delay_pmf <- function(neg_binom_mu = 6.98665,
                                          neg_binom_size = 2.490848) {
  density <- dnbinom(
    x = seq(0, 30, 1),
    mu = neg_binom_mu, size = neg_binom_size
  )
  hosp_onset_delay_pmf <- density / sum(density)

  return(hosp_onset_delay_pmf)
}


#' @title Make reporting delay pmf
#' @description
#' Convolve the incubation period pmf with the symptom to hospital admission pmf
#' and normalize
#'
#' @param incubation_period_pmf a numeric vector, sums to 1, indicating
#' the probability of time from infection to symptom onset
#' @param hospital_onset_delay_pmf a numeric vector, sums to 1, indicating the
#' proabbility of time from symptom onset to hospital admissions
#'
#' @return convolution of incubation period and sympton onset to hospital
#' admission pmf
make_reporting_delay_pmf <- function(incubation_period_pmf,
                                     hospital_onset_delay_pmf) {
  pmfs <- list(
    "incubation_period" = incubation_period_pmf,
    "hosp_onset_delay" = hospital_onset_delay_pmf
  )

  infection_to_hosp_delay_pmf <- add_pmfs(pmfs) |>
    wwinference::to_simplex()
  return(infection_to_hosp_delay_pmf)
}
