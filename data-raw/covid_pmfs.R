# Load in the parameters
params <- wwinference::get_params(
  system.file("extdata", "example_params.toml",
    package = "wwinference"
  )
)

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
