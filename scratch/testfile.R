library(wwinference)
library(dplyr)
library(ggplot2)
library(tidybayes)
library(tidyverse)

hosp_data <- wwinference::hosp_data
hosp_data_eval <- wwinference::hosp_data_eval
ww_data <- wwinference::ww_data
rt_global <- wwinference::rt_global_data
rt_site <- wwinference::rt_site_data

head(ww_data)
head(hosp_data)


params <- get_params(
  system.file("extdata", "example_params.toml",
    package = "wwinference"
  )
)


ww_data_preprocessed <- wwinference::preprocess_ww_data(
  ww_data,
  conc_col_name = "genome_copies_per_ml",
  lod_col_name = "lod"
)


hosp_data_preprocessed <- wwinference::preprocess_hosp_data(
  hosp_data,
  count_col_name = "daily_hosp_admits",
  pop_size_col_name = "state_pop"
)

# Sites ------------------------------------------------------------------------
site_locs <- data.frame(
  x = c(85, 37, 36, 7),
  y = c(12, 75, 75, 96),
  ID = c("Site 1", "Site 2", "Site 3", "Site 4")
)
ggplot(data = site_locs) +
  geom_point(aes(x = x, y = y, colour = ID, shape = ID), size = 5) +
  labs(title = "Fake Locations") +
  guides(
    colour = guide_legend(title = "Locations"),
    shape = guide_legend(title = "Locations")
  ) +
  theme_bw()

#-------------------------------------------------------------------------------
ggplot(ww_data_preprocessed) +
  geom_point(
    aes(
      x = date, y = genome_copies_per_ml,
      color = as.factor(lab_site_name)
    ),
    show.legend = FALSE
  ) +
  geom_point(
    data = ww_data_preprocessed |> filter(genome_copies_per_ml <= lod),
    aes(x = date, y = genome_copies_per_ml, color = "red"),
    show.legend = FALSE
  ) +
  geom_hline(aes(yintercept = lod), linetype = "dashed") +
  facet_wrap(~lab_site_name, scales = "free") +
  xlab("") +
  ylab("Genome copies/mL") +
  ggtitle("Lab-site level wastewater concentration") +
  theme_bw()

ggplot(hosp_data_preprocessed) +
  # Plot the hospital admissions data that we will evaluate against in white
  geom_point(
    data = hosp_data_eval, aes(
      x = date,
      y = daily_hosp_admits_for_eval
    ),
    shape = 21, color = "black", fill = "white"
  ) +
  # Plot the data we will calibrate to
  geom_point(aes(x = date, y = count)) +
  xlab("") +
  ylab("Daily hospital admissions") +
  ggtitle("Global level hospital admissions") +
  theme_bw()

# Rt site-----------------------------------------------------------------------
rt_site_df <- as.data.frame(t(rt_site)) %>%
  mutate(date = hosp_data_eval$date) %>%
  `colnames<-`(c(
    "Site: 1",
    "Site: 2",
    "Site: 3",
    "Site: 4",
    "remainder of pop",
    "date"
  )) %>%
  pivot_longer(cols = -date) %>%
  `colnames<-`(c(
    "date",
    "subpop",
    "value"
  ))
ggplot(rt_site_df) +
  geom_line(aes(
    x = date,
    y = value
  )) +
  facet_wrap(~subpop, scales = "free") +
  xlab("") +
  ylab("Subpopulation R(t)") +
  ggtitle("R(t) estimate") +
  theme_bw()
#-------------------------------------------------------------------------------
# Rt global---------------------------------------------------------------------
rt_global_df <- as.data.frame(rt_global) %>%
  mutate(date = hosp_data_eval$date)
ggplot(rt_global_df) +
  geom_line(aes(
    x = date,
    y = rt_global
  )) +
  xlab("") +
  ylab("Global R(t)") +
  ggtitle("Global R(t) estimate") +
  theme_bw()
#-------------------------------------------------------------------------------

ww_data_to_fit <- wwinference::indicate_ww_exclusions(
  ww_data_preprocessed,
  outlier_col_name = "flag_as_ww_outlier",
  remove_outliers = TRUE
)


forecast_date <- "2023-12-06"
calibration_time <- 90
forecast_horizon <- 28


generation_interval <- wwinference::generation_interval
inf_to_hosp <- wwinference::inf_to_hosp

# Assign infection feedback equal to the generation interval
infection_feedback_pmf <- generation_interval


model <- wwinference::compile_model(
  model_filepath = "inst/stan/wwinference.stan",
  include_paths = "inst/stan"
)


fit <- wwinference::wwinference(
  ww_data_to_fit,
  hosp_data_preprocessed,
  model_spec = get_model_spec(
    forecast_date = forecast_date,
    calibration_time = calibration_time,
    forecast_horizon = forecast_horizon,
    generation_interval = generation_interval,
    inf_to_count_delay = inf_to_hosp,
    infection_feedback_pmf = infection_feedback_pmf
  ),
  mcmc_options = get_mcmc_options(),
  compiled_model = model
)


head(fit)


draws_df <- fit$draws_df
sampled_draws <- sample(1:max(draws_df$draw), 100)

# Hospital admissions: fits, nowcasts, forecasts
ggplot(draws_df |> dplyr::filter(
  name == "pred_counts",
  draw %in% sampled_draws
)) +
  geom_line(aes(x = date, y = pred_value, group = draw),
    color = "red4", alpha = 0.1, size = 0.2
  ) +
  geom_point(
    data = hosp_data_eval,
    aes(x = date, y = daily_hosp_admits_for_eval),
    shape = 21, color = "black", fill = "white"
  ) +
  geom_point(aes(x = date, y = observed_value)) +
  geom_vline(aes(xintercept = lubridate::ymd(forecast_date)),
    linetype = "dashed"
  ) +
  xlab("") +
  ylab("Daily hospital admissions") +
  ggtitle("Fit and forecasted hospital admissions ") +
  theme_bw()

# R(t) of the hypothetical state
ggplot(draws_df |> dplyr::filter(
  name == "global R(t)",
  draw %in% sampled_draws
)) +
  geom_line(data = rt_global_df, aes(
    x = rt_global_df$date,
    y = rt_global_df$rt_global
  ), color = "red") +
  ggplot2::geom_step(
    aes(x = date, y = pred_value, group = draw),
    color = "blue4",
    alpha = 0.1, linewidth = 0.2
  ) +
  # geom_line(aes(x = date, y = pred_value, group = draw),
  #  color = "blue4", alpha = 0.1, size = 0.2
  # ) +
  geom_vline(aes(xintercept = lubridate::ymd(forecast_date)),
    linetype = "dashed"
  ) +
  geom_hline(aes(yintercept = 1), linetype = "dashed") +
  xlab("") +
  ylab("Global R(t) ") +
  ggtitle("Global R(t) estimate (Red line is actual)") +
  theme_bw()


ggplot(draws_df |> dplyr::filter(
  name == "pred_ww",
  draw %in% sampled_draws
) |>
  dplyr::mutate(
    site_lab_name = glue::glue("{subpop}, Lab: {lab}")
  )) +
  geom_line(
    aes(
      x = date, y = log(pred_value),
      color = subpop,
      group = draw
    ),
    alpha = 0.1, size = 0.2,
    show.legend = FALSE
  ) +
  geom_point(aes(x = date, y = log(observed_value)),
    color = "black", show.legend = FALSE
  ) +
  facet_wrap(~site_lab_name, scales = "free") +
  geom_vline(aes(xintercept = lubridate::ymd(forecast_date)),
    linetype = "dashed"
  ) +
  xlab("") +
  ylab("Log(Genome copies/mL)") +
  ggtitle("Lab-site level wastewater concentration") +
  theme_bw()

ggplot(draws_df |> dplyr::filter(
  name == "subpop R(t)",
  draw %in% sampled_draws
)) +
  geom_line(data = rt_site_df, aes(
    x = date,
    y = value
  ), color = "red") +
  # //geom_line(
  # // aes(
  # //   x = date, y = pred_value, group = draw,
  # //   color = subpop
  # // ),
  # // alpha = 0.1, size = 0.2
  # //) +
  ggplot2::geom_step(
    aes(x = date, y = pred_value, group = draw, color = subpop),
    alpha = 0.1, linewidth = 0.2
  ) +
  geom_vline(aes(xintercept = lubridate::ymd(forecast_date)),
    linetype = "dashed"
  ) +
  facet_wrap(~subpop, scales = "free") +
  geom_hline(aes(yintercept = 1), linetype = "dashed") +
  xlab("") +
  ylab("Subpopulation R(t)") +
  ggtitle("Site R(t) estimate (Red line is actual)") +
  theme_bw()



summary_spatial_params <- fit$raw_fit_obj$summary() %>%
  filter(variable %in% c(
    "autoreg_rt_site",
    "phi",
    "sigma_generalized",
    "scaling_factor"
  )) %>%
  mutate(actual_values = c(0.6, 0.2, 0.05^4, 1.1))
summary_spatial_params




temp_draws <- fit$raw_fit_obj$draws(variables = c(
  "autoreg_rt_site",
  "phi",
  "sigma_generalized",
  "scaling_factor"
))
temp_draws_df <- as.data.frame(as.table(temp_draws))
names(temp_draws_df) <- c("iteration", "chain", "variable", "value")
temp_draws_df <- temp_draws_df %>%
  mutate(
    variable = case_when(
      variable == "autoreg_rt_site" ~ "AR Coefficient on Delta Terms",
      variable == "phi" ~ "Phi for Exp. Corr. Func.",
      variable == "sigma_generalized" ~ "Generalized Variance",
      variable == "scaling_factor" ~ "ScalingFactor"
    ),
    variable = factor(variable),
  )
actual_values <- c(
  "AR Coefficient on Delta Terms" = 0.6,
  "Phi for Exp. Corr. Func." = 0.2,
  "Generalized Variance" = 0.00000625,
  "ScalingFactor" = 1.1
)
actual_values_df <- data.frame(
  variable = names(actual_values),
  actual_value = as.vector(actual_values)
)
temp_draws_df <- temp_draws_df %>%
  left_join(
    actual_values_df,
    by = "variable"
  )
ggplot(
  data = temp_draws_df
) +
  geom_histogram(
    aes(
      x = value
    ),
    color = "white",
    fill = "darkblue"
  ) +
  geom_vline(
    aes(xintercept = actual_value),
    color = "red2",
    linetype = "dashed",
    size = 1.5
  ) +
  facet_grid(~variable, scales = "free") +
  xlab("Sampled Value") +
  ylab("count") +
  ggtitle(
    "Histograms of Spatial Parameters (red line is actual simulation value)"
  ) +
  theme( # //axis.title.x = element_blank(),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 14),
    axis.line = element_line(colour = "black", size = 1.25),
    axis.ticks = element_line(colour = "black", size = 1.5),
    axis.ticks.length = unit(.25, "cm"),
    panel.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(colour = "white", face = "bold", size = 12),
    legend.background = element_rect(fill = "black"),
    legend.key.width = unit(.025, "npc"),
    plot.title = element_text(face = "bold", size = 16),
    strip.text = element_text(colour = "white", face = "bold", size = 12),
    strip.background = element_rect(fill = "black")
  )




plot_hosp <- get_plot_forecasted_counts(
  draws = draws_df,
  count_data_eval = hosp_data_eval,
  count_data_eval_col_name = "daily_hosp_admits_for_eval",
  forecast_date = forecast_date
)
plot_hosp
plot_ww <- get_plot_ww_conc(draws_df, forecast_date)
plot_ww
plot_state_rt <- get_plot_global_rt(draws_df, forecast_date)
plot_state_rt
plot_subpop_rt <- get_plot_subpop_rt(draws_df, forecast_date)
plot_subpop_rt


convergence_flag_df <- wwinference::get_model_diagnostic_flags(
  stan_fit_obj =
    fit$raw_fit_obj
)
parameter_diangostics <- fit$raw_fit_obj$summary()
