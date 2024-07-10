# Load in the parameters
params <- get_params(
  system.file("extdata", "example_params.toml",
    package = "wwinference"
  )
)




# Simulate daily double censored PMF. From {epinowcast}:
# https://package.epinowcast.org/dev/reference/simulate_double_censored_pmf.html #nolint
#
# This function simulates the probability mass function of a  daily
# double-censored process. The process involves two distributions: a primary
# distribution which represents the censoring process for the primary event
# and another distribution (which is offset by the primary).
#
# Based off of:
# https://www.medrxiv.org/content/10.1101/2024.01.12.24301247v1
simulate_double_censored_pmf <- function(
    max, fun_primary = stats::runif, primary_args = list(),
    fun_dist = stats::rlnorm,
    dist_args = list(...), n = 1e6, ...) {
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


# This retruns an incubation period pmf corresponding to the incubation period
# for COVID after Omicron used in
# Park et al 2023. These estimates are from early Omicron.

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

  inc_period_pmf <- to_simplex(discr_gr_adj_weibull$density0)
  return(inc_period_pmf)
}

# Makes the hospital onset delay pmf using the parameter estimates based on
# Danache et al linelist data from symptom onset to hospital
# admission. See below:
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0261428

make_hospital_onset_delay_pmf <- function(neg_binom_mu = 6.98665,
                                          neg_binom_size = 2.490848) {
  density <- dnbinom(
    x = seq(0, 30, 1),
    mu = neg_binom_mu, size = neg_binom_size
  )
  hosp_onset_delay_pmf <- density / sum(density)

  return(hosp_onset_delay_pmf)
}


# Makes the reporting delay pmf by convolving the incubation period pmf with
# the symptom to hospital admission pmf and normalizing

make_reporting_delay_pmf <- function(incubation_period_pmf,
                                     hospital_onset_delay_pmf) {
  pmfs <- list(
    "incubation_period" = incubation_period_pmf,
    "hosp_onset_delay" = hospital_onset_delay_pmf
  )

  infection_to_hosp_delay_pmf <- add_pmfs(pmfs) |>
    to_simplex()
  return(infection_to_hosp_delay_pmf)
}





# Put it all together
generation_interval <- withr::with_seed(42, {
  wwinference::simulate_double_censored_pmf(
    max = params$gt_max, meanlog = params$mu_gi,
    sdlog = params$sigma_gi, fun_dist = rlnorm, n = 5e6
  ) |> wwinference::drop_first_and_renormalize()
})

inc <- wwinference::make_incubation_period_pmf(
  params$backward_scale, params$backward_shape, params$r
)
sym_to_hosp <- wwinference::make_hospital_onset_delay_pmf(
  params$neg_binom_mu,
  params$neg_binom_size
)
inf_to_hosp <- wwinference::make_reporting_delay_pmf(inc, sym_to_hosp)

usethis::use_data(generation_interval, overwrite = TRUE)
usethis::use_data(inf_to_hosp, overwrite = TRUE)
