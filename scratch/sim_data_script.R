# Specify all of the simulation settings (eventually we will put this into
# vignette_data.R but leave as scratch for now)

# Hacky way of pulling in all the functions in R
list.files(file.path("R"), full.names = TRUE) |>
  purrr::walk(source)

r_in_weeks <- c(
  rep(1.1, 5), rep(0.9, 5),
  1 + 0.007 * 1:16
)
n_sites <- 4
ww_pop_sites <- c(4e5, 2e5, 1e5, 5e4)
pop_size <- 3e6
site <- c(1, 1, 2, 3, 4)
lab <- c(1, 2, 3, 3, 3)
ot <- 90
nt <- 9
forecast_horizon <- 28
sim_start_date <- lubridate::ymd("2023-09-01")
hosp_wday_effect <- c(
  0.95, 1.01, 1.02,
  1.02, 1.01, 1,
  0.99
) / 7
i0_over_n <- 5e-4
initial_growth <- 1e-4
sd_in_lab_level_multiplier <- 0.25
mean_obs_error_in_ww_lab_site <- 0.2
mean_reporting_freq <- 1 / 5
sd_reporting_freq <- 1 / 20
mean_reporting_latency <- 7
sd_reporting_latency <- 3
mean_log_lod <- 3.8
sd_log_lod <- 0.2
global_rt_sd <- 0.03
sigma_eps <- 0.05
sd_i0_over_n <- 0.5
infection_feedback <- TRUE
subpop_phi <- c(25, 50, 70, 40, 100)
input_params_path <-
  fs::path_package("extdata",
    "example_params.toml",
    package = "wwinference"
  )
if_feedback <- FALSE

# R(t) comparison
rt_r <- new_i_over_n / (convolve(new_i_over_n,
  rev(c(0, generation_interval)),
  type = "open"
))[1:(uot + ot + ht)]

test <- tibble::tibble(
  rt_stan = rt,
  rt_r = rt_r,
  t = 1:(ot + ht)
)

ggplot(test) +
  geom_line(aes(x = t, y = rt_r),
    color = "black",
    linewidth = 2
  ) +
  geom_line(aes(x = t, y = rt_stan), color = "red")


new_i_test <- rt_r * (convolve(new_i_over_n, rev(generation_interval), type = "open")[1:(ot + ht)]) # nolint
plot(new_i_test)
