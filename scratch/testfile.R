# Expose the stan functions into the global environment
model <- cmdstanr::cmdstan_model(
  stan_file = file.path("inst", "stan", "wwinference.stan"),
  compile = TRUE,
  compile_standalone = TRUE,
  force_recompile = TRUE
)
model$expose_functions(global = TRUE)
model <- cmdstanr::cmdstan_model(
  stan_file = file.path("inst", "stan", "functions", "spatial_functions.stan"),
  compile = TRUE,
  compile_standalone = TRUE,
  force_recompile = TRUE
)
model$expose_functions(global = TRUE)

n_time <- 150
state_dev_ar_coeff <- 0.8
log_state_rt <- rnorm(
  n = n_time,
  mean = 1.2,
  sd = 0.05
)
state_dev_noise_vec <- state_deviation_noise_vec_aux_rng(
  scaling_factor = 1.1,
  sigma_eps = sqrt(0.2),
  n_time = n_time
)
stan_log_aux_site_rt <- construct_aux_rt(
  log_state_rt = log_state_rt,
  state_deviation_ar_coeff = state_dev_ar_coeff,
  state_deviation_noise_vec = state_dev_noise_vec
)


state_deviation_t_i <- 0
log_aux_site_rt <- matrix(
  data = 0,
  ncol = n_time,
  nrow = 1
)
for (t_i in 1:n_time) {
  state_deviation_t_i <- state_dev_ar_coeff * state_deviation_t_i +
    state_dev_noise_vec[t_i]
  log_aux_site_rt[t_i] <- log_state_rt[t_i] + state_deviation_t_i
}

stan_log_aux_site_rt == log_aux_site_rt
