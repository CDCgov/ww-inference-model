# This script shows how to generate each of the plots from the Model Run
# Studio output

library(ggplot2)
library(readr)
library(dplyr)
library(lubridate)
library(scales)

# Set-up--------------------------------------------------------------------
n_draws_to_plot <- 100

# Load in the posterior draws dataframes-------------------------------------
predicted_hosp <- read_csv("predicted_hosp.csv")
predicted_ww <- read_csv("predicted_ww.csv")
global_rt <- read_csv("global_rt.csv")
subpop_rt <- read_csv("subpop_rt.csv")

sampled_draws <- sample(1:max(predicted_hosp$draw), n_draws_to_plot)

# Predicted hospital admissions--------------------------------------------
hosp_draws_to_plot <- predicted_hosp |>
  filter(draw %in% !!sampled_draws)

plot_hosp_forecast <- ggplot(draws_to_plot) +
  geom_line(
    aes(x = date, y = pred_value, group = draw),
    color = "red4", alpha = 0.1, linewidth = 0.2
  ) +
  geom_vline(
    xintercept = ymd(forecast_date),
    linetype = "dashed"
  ) +
  geom_point(aes(x = date, y = observed_value)) +
  xlab("") +
  ylab("Daily hospital admissions") +
  ggtitle("Fit and forecasted hospital admissions") +
  scale_x_date(
    date_breaks = "2 weeks",
    labels = date_format("%Y-%m-%d")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      size = 8, vjust = 1,
      hjust = 1, angle = 45
    ),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(
      size = 10,
      vjust = 0.5, hjust = 0.5
    )
  )

plot_hosp_forecast

# Predicted wastewater concentrations --------------------------------------

ww_draws_to_plot <- predicted_ww |>
  filter(draw %in% !!sampled_draws)

plot_ww_forecast <- ggplot(ww_draws_to_plot) +
  geom_line(
    aes(
      x = date, y = pred_value,
      color = subpop_name,
      group = draw
    ),
    alpha = 0.1, size = 0.2,
    show.legend = FALSE
  ) +
  geom_point(aes(x = date, y = observed_value),
    color = "black", show.legend = FALSE, size = 0.5
  ) +
  geom_point(
    data = ww_draws_to_plot |> filter(below_lod == 1),
    aes(x = date, y = observed_value),
    color = "blue", show.legend = FALSE, size = 0.5
  ) +
  facet_wrap(~lab_site_name, scales = "free_y") +
  geom_vline(
    xintercept = ymd(forecast_date),
    linetype = "dashed"
  ) +
  xlab("") +
  ylab("Log genome copies/mL") +
  ggtitle("Lab-site level wastewater concentrations") +
  scale_x_date(
    date_breaks = "2 weeks",
    labels = date_format("%Y-%m-%d")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      size = 5, vjust = 1,
      hjust = 1, angle = 45
    ),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 5),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 6),
    plot.title = element_text(
      size = 10,
      vjust = 0.5, hjust = 0.5
    )
  )

plot_ww_forecast

# Global R(t) ------------------------------------------------------------
global_rt_draws_to_plot <- global_rt |>
  filter(draw %in% !!sampled_draws)

plot_global_rt <- ggplot(global_rt_draws_to_plot) +
  geom_step(
    aes(x = date, y = pred_value, group = draw),
    color = "blue4", alpha = 0.1, linewidth = 0.2
  ) +
  geom_vline(
    xintercept = ymd(forecast_date),
    linetype = "dashed"
  ) +
  geom_hline(aes(yintercept = 1), linetype = "dashed") +
  xlab("") +
  ylab("Global R(t)") +
  ggtitle("Global R(t) estimate") +
  scale_x_date(
    date_breaks = "2 weeks",
    labels = date_format("%Y-%m-%d")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      size = 5, vjust = 1,
      hjust = 1, angle = 45
    ),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 5),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(
      size = 10,
      vjust = 0.5, hjust = 0.5
    )
  )

plot_global_rt

# Subpopulation R(t) ---------------------------------------------------------
subpop_rt_draws_to_plot <- subpop_rt |>
  filter(draw %in% !!sampled_draws)

plot_subpop_rt <- ggplot(subpop_rt_draws_to_plot) +
  geom_step(
    aes(
      x = date, y = pred_value, group = draw,
      color = subpop_name
    ),
    alpha = 0.1, linewidth = 0.2,
    show.legend = FALSE
  ) +
  geom_vline(
    xintercept = ymd(forecast_date),
    linetype = "dashed",
    show.legend = FALSE
  ) +
  facet_wrap(~subpop_name) +
  geom_hline(aes(yintercept = 1), linetype = "dashed") +
  xlab("") +
  ylab("Subpopulation R(t)") +
  ggtitle("R(t) estimate of each subpopulation") +
  scale_x_date(
    date_breaks = "2 weeks",
    labels = date_format("%Y-%m-%d")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      size = 5, vjust = 1,
      hjust = 1, angle = 45
    ),
    axis.text.y = element_text(size = 5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 6),
    plot.title = element_text(
      size = 10,
      vjust = 0.5, hjust = 0.5
    )
  )

plot_subpop_rt
