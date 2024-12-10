## Investigating if prior is identifiable / effect of changing prior
# Library
library(wwinference)
library(tidyverse)
library(patchwork)


## Presets - Data Generation
set.seed(2024)
n_sites <- 4
fake_locs <- data.frame(
  x = runif(n_sites, 0, 100),
  y = runif(n_sites, 0, 100),
  ID = paste0("Site ", 1:n_sites)
)
r_in_weeks <- c(
  rep(1.1, 5), rep(0.9, 5),
  0.9 + 0.0005 * (1:16)^2
)
site <- 1:n_sites
lab <- rep(1, n_sites)
ww_pop_sites <- c(5e5, 7.5e5, 3.8e5, 4.5e5)
pop_size <- 2 * sum(ww_pop_sites)
global_rt_sd <- 0
exp_corr_param <- list(
  dist_matrix = as.matrix(
    dist(
      data.frame(
        x = fake_locs$x,
        y = fake_locs$y
      ),
      diag = TRUE,
      upper = TRUE
    )
  ),
  phi = 1.5,
  l = 1
)
sigma_sqrd_generalized <- 0.01^n_sites

# Data Generation - Exponential
simulated_data <- withr::with_seed(1, {
  wwinference::generate_simulated_data(
    r_in_weeks = r_in_weeks, # nolint
    n_sites = n_sites,
    site = site,
    lab = lab,
    ww_pop_sites = ww_pop_sites,
    pop_size = pop_size,
    global_rt_sd = global_rt_sd,
    use_spatial_corr = TRUE,
    corr_function = exponential_decay_corr_func_r,
    corr_fun_params = exp_corr_param,
    sigma_sqrd_generalized = sigma_sqrd_generalized,
    aux_site_bool = TRUE
  )
})
hosp_data_preprocessed <- wwinference::preprocess_count_data(
  simulated_data$hosp_data,
  count_col_name = "daily_hosp_admits",
  pop_size_col_name = "state_pop"
)
ww_data_preprocessed <- wwinference::preprocess_ww_data(
  simulated_data$ww_data
)
ww_data_to_fit <- wwinference::indicate_ww_exclusions(
  ww_data_preprocessed,
  outlier_col_name = "flag_as_ww_outlier",
  remove_outliers = TRUE
)


## Presets - Fit
params_new <- get_params(
  system.file("extdata", "example_params.toml",
    package = "wwinference"
  )
)
cat("New Log Phi Mean :", params_new$log_phi_mu_prior, "\n")
params_old <- params_new
params_old$log_phi_mu_prior <- log(0.25)
params_old$log_phi_sd_prior <- 0.75
cat("Old Log Phi Mean :", params_old$log_phi_mu_prior, "\n")
forecast_date <- "2023-12-06"
calibration_time <- 90
forecast_horizon <- 28

generation_interval <- wwinference::default_covid_gi
inf_to_hosp <- wwinference::default_covid_inf_to_hosp

# Assign infection feedback equal to the generation interval
infection_feedback_pmf <- generation_interval

model <- wwinference::compile_model()

iter_warmup <- 250
iter_sampling <- 500
fit_options <- get_mcmc_options(
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling,
  seed = 123
)


## Fitting
fit_old <- wwinference::wwinference(
  ww_data = ww_data_to_fit,
  count_data = hosp_data_preprocessed,
  forecast_date = forecast_date,
  calibration_time = calibration_time,
  forecast_horizon = forecast_horizon,
  model_spec = get_model_spec(
    generation_interval = generation_interval,
    inf_to_count_delay = inf_to_hosp,
    infection_feedback_pmf = infection_feedback_pmf,
    params = params_old
  ),
  fit_opts = fit_options,
  compiled_model = model,
  dist_matrix = as.matrix(exp_corr_param$dist_matrix),
  corr_structure_switch = 1
)
phi_draws_old <- fit_old$fit$result$draws(
  variables = "phi",
  format = "draws_array"
)
phi_draws_old <- phi_draws_old %>%
  as.data.frame() %>%
  pivot_longer(
    everything(),
    names_to = "chain",
    values_to = "phi"
  ) %>%
  mutate(
    chain = case_when(
      chain == "1.phi" ~ 1,
      chain == "2.phi" ~ 2,
      chain == "3.phi" ~ 3,
      chain == "4.phi" ~ 4
    ),
    chain = factor(chain)
  )

fit_new <- wwinference::wwinference(
  ww_data = ww_data_to_fit,
  count_data = hosp_data_preprocessed,
  forecast_date = forecast_date,
  calibration_time = calibration_time,
  forecast_horizon = forecast_horizon,
  model_spec = get_model_spec(
    generation_interval = generation_interval,
    inf_to_count_delay = inf_to_hosp,
    infection_feedback_pmf = infection_feedback_pmf,
    params = params_new
  ),
  fit_opts = fit_options,
  compiled_model = model,
  dist_matrix = as.matrix(exp_corr_param$dist_matrix),
  corr_structure_switch = 1
)
phi_draws_new <- fit_new$fit$result$draws(
  variables = "phi",
  format = "draws_array"
)
phi_draws_new <- phi_draws_new %>%
  as.data.frame() %>%
  pivot_longer(
    everything(),
    names_to = "chain",
    values_to = "phi"
  ) %>%
  mutate(
    chain = case_when(
      chain == "1.phi" ~ 1,
      chain == "2.phi" ~ 2,
      chain == "3.phi" ~ 3,
      chain == "4.phi" ~ 4
    ),
    chain = factor(chain)
  )



## Visualizing Fits
# Trace Plots
phi_trace_new_plot <- ggplot(phi_draws_new) +
  geom_line(
    aes(
      y = phi,
      x = seq_along(phi),
      colour = chain
    )
  ) +
  labs(
    title = "Trace Plot for Phi Draws w/ New Prior",
    x = "iter",
    y = "phi",
    colour = "Chain"
  ) +
  scale_color_manual(
    values = c("gold", "darkred", "darkgreen", "darkblue")
  ) +
  theme_bw()
phi_trace_old_plot <- ggplot(phi_draws_old) +
  geom_line(
    aes(
      y = phi,
      x = seq_along(phi),
      colour = chain
    )
  ) +
  labs(
    title = "Trace Plot for Phi Draws w/ Old Prior",
    x = "iter",
    y = "phi",
    colour = "Chain"
  ) +
  scale_color_manual(
    values = c("gold", "darkred", "darkgreen", "darkblue")
  ) +
  theme_bw()

phi_trace_old_plot / phi_trace_new_plot +
  plot_layout(
    guides = "collect"
  )

# Prior and Posterior
prior_old <- exp(rnorm(
  n = 1000,
  mean = params_old$log_phi_mu_prior,
  sd = params_old$log_phi_sd_prior
))
prior_new <- exp(rnorm(
  n = 1000,
  mean = params_new$log_phi_mu_prior,
  sd = params_new$log_phi_sd_prior
))

prior_post_old_plot <- ggplot() +
  geom_vline(
    xintercept = exp_corr_param$phi,
    linetype = "dashed"
  ) +
  geom_density(
    aes(
      x = prior_old,
      fill = "Prior",
      colour = "Prior"
    ),
    alpha = .25
  ) +
  geom_density(
    data = phi_draws_old %>%
      filter(
        chain == 1
      ),
    aes(
      x = phi,
      fill = "Posterior",
      colour = "Posterior"
    ),
    alpha = .25
  ) +
  scale_color_manual(
    breaks = c("Prior", "Posterior"),
    values = c("dodgerblue", "purple")
  ) +
  scale_fill_manual(
    breaks = c("Prior", "Posterior"),
    values = c("dodgerblue", "purple")
  ) +
  labs(
    title = "Prior & Posterior w/ Old Prior",
    subtitle = "Vertical Dashed Line Actual",
    fill = "Dist. Type",
    colour = "Dist. Type",
    y = "Density",
    x = "Value"
  ) +
  theme_bw()

prior_post_new_plot <- ggplot() +
  geom_vline(
    xintercept = exp_corr_param$phi,
    linetype = "dashed"
  ) +
  geom_density(
    aes(
      x = prior_new,
      fill = "Prior",
      colour = "Prior"
    ),
    alpha = .25
  ) +
  geom_density(
    data = phi_draws_new %>%
      filter(
        chain == 1
      ),
    aes(
      x = phi,
      fill = "Posterior",
      colour = "Posterior"
    ),
    alpha = .25
  ) +
  scale_color_manual(
    breaks = c("Prior", "Posterior"),
    values = c("dodgerblue", "purple")
  ) +
  scale_fill_manual(
    breaks = c("Prior", "Posterior"),
    values = c("dodgerblue", "purple")
  ) +
  labs(
    title = "Prior & Posterior w/ New Prior",
    subtitle = "Vertical Dashed Line Actual",
    fill = "Dist. Type",
    colour = "Dist. Type",
    y = "Density",
    x = "Value"
  ) +
  theme_bw()

prior_post_old_plot / prior_post_new_plot +
  plot_layout(
    guides = "collect"
  )
