library(wwinference)
library(dplyr)
library(ggplot2)
library(tidybayes)


hosp_data <- wwinference::hosp_data
hosp_data_eval <- wwinference::hosp_data_eval
ww_data <- wwinference::ww_data

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
  ggtitle("State level hospital admissions") +
  theme_bw()


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
  geom_line(aes(x = date, y = pred_value, group = draw),
    color = "blue4", alpha = 0.1, size = 0.2
  ) +
  geom_vline(aes(xintercept = lubridate::ymd(forecast_date)),
    linetype = "dashed"
  ) +
  geom_hline(aes(yintercept = 1), linetype = "dashed") +
  xlab("") +
  ylab("Global R(t) ") +
  ggtitle("Global R(t) estimate") +
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
  geom_line(
    aes(
      x = date, y = pred_value, group = draw,
      color = subpop
    ),
    alpha = 0.1, size = 0.2
  ) +
  geom_vline(aes(xintercept = lubridate::ymd(forecast_date)),
    linetype = "dashed"
  ) +
  facet_wrap(~subpop, scales = "free") +
  geom_hline(aes(yintercept = 1), linetype = "dashed") +
  xlab("") +
  ylab("Subpopulation R(t)") +
  ggtitle("R(t) estimate") +
  theme_bw()


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
