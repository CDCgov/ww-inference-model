# testfile
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



hosp_data_preprocessed <- wwinference::preprocess_count_data(
  hosp_data,
  count_col_name = "daily_hosp_admits",
  pop_size_col_name = "state_pop"
)


ww_data_preprocessed <- wwinference::preprocess_ww_data(
  ww_data,
  conc_col_name = "genome_copies_per_ml",
  lod_col_name = "lod"
)

# Sites ------------------------------------------------------------------------
site_locs <- data.frame(
  x = c(85, 37, 36, 7),
  y = c(12, 75, 75, 96),
  ID = c("Site 1", "Site 2", "Site 3", "Site 4")
)
ggplot(data = site_locs) +
  geom_point(aes(x = x, y = y, colour = ID, shape = ID), size = 15) +
  labs(title = "Fake Locations") +
  guides(
    colour = guide_legend(title = "Locations"),
    shape = guide_legend(title = "Locations")
  ) +
  theme_bw()
dist_matrix <- dist(
  site_locs[, c(1, 2)],
  upper = TRUE,
  diag = TRUE
)

# Correlation matrix -----------------------------------------------------------
corr_function <- exponential_decay_corr_func_r
corr_fun_params <- list(
  dist_matrix = as.matrix(
    dist(
      data.frame(
        x = c(85, 37, 36, 7),
        y = c(12, 75, 75, 96)
      ),
      diag = TRUE,
      upper = TRUE
    )
  ) / 114.62984,
  phi = 0.2,
  l = 1
)
corr_matrix <- corr_function(corr_fun_params)
corr_df <- as.data.frame(corr_matrix) %>%
  `colnames<-`(c(
    "Site 1",
    "Site 2",
    "Site 3",
    "Site 4"
  )) %>%
  `rownames<-`(c(
    "Site 1",
    "Site 2",
    "Site 3",
    "Site 4"
  )) %>%
  rownames_to_column(var = "Var1") %>%
  pivot_longer(
    -Var1,
    names_to = "Var2",
    values_to = "value"
  )
ggplot(corr_df, aes(Var1, Var2)) +
  geom_point(
    aes(size = abs(value), fill = value),
    shape = 21,
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1),
    space = "Lab",
    name = "Correlation"
  ) +
  scale_size_continuous(
    range = c(1, 20),
    guide = "none"
  ) +
  coord_fixed() +
  ylab("") +
  xlab("") +
  ggtitle("Exponential Correlation Matrix Visual") +
  theme( # //axis.title.x = element_blank(),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 14),
    axis.line = element_line(colour = "black", size = 1.25),
    axis.ticks = element_line(colour = "black", size = 1.5),
    axis.ticks.length = unit(.25, "cm"),
    panel.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(colour = "white", face = "bold", size = 12),
    legend.text = element_text(colour = "white", face = "bold", size = 12),
    legend.background = element_rect(fill = "black"),
    legend.key.width = unit(.025, "npc"),
    plot.title = element_text(face = "bold", size = 16),
    strip.text = element_text(colour = "white", face = "bold", size = 12),
    strip.background = element_rect(fill = "black")
  )
#-------------------------------------------------------------------------------

ww_data_to_fit <- wwinference::indicate_ww_exclusions(
  ww_data_preprocessed,
  outlier_col_name = "flag_as_ww_outlier",
  remove_outliers = TRUE
)


forecast_date <- "2023-12-06"
calibration_time <- 90
forecast_horizon <- 28


generation_interval <- wwinference::default_covid_gi
inf_to_hosp <- wwinference::default_covid_inf_to_hosp

# Assign infection feedback equal to the generation interval
infection_feedback_pmf <- generation_interval


model <- wwinference::compile_model(
  model_filepath = "inst/stan/wwinference.stan",
  include_paths = "inst/stan"
)


fit <- wwinference::wwinference(
  ww_data = ww_data_to_fit,
  count_data = hosp_data_preprocessed,
  forecast_date = forecast_date,
  calibration_time = calibration_time,
  forecast_horizon = forecast_horizon,
  model_spec = get_model_spec(
    generation_interval = generation_interval,
    inf_to_count_delay = inf_to_hosp,
    infection_feedback_pmf = infection_feedback_pmf,
    params = params
  ),
  mcmc_options = get_mcmc_options(),
  compiled_model = model,
  dist_matrix = as.matrix(dist_matrix),
  corr_func = "ljk"
  # // dist_matrix = NULL,
  # // bool_spatial_comp = TRUE
)


head(fit)


draws_df <- fit$draws_df
sampled_draws <- sample(1:max(draws_df$draw), 100)


# Correlation medians from stan output
cor_matrix_draws <- fit$raw_fit_obj$draws(
  variables = "non_norm_omega",
  format = "draws_array"
)
median_cor_matrix <- apply(
  cor_matrix_draws,
  c(2, 3),
  median
)
median_cor_matrix_df <- as.table(median_cor_matrix) %>%
  as.data.frame() %>%
  filter(
    chain == 1
  ) %>%
  select(
    -chain
  )
median_cor_matrix_df <- median_cor_matrix_df %>%
  mutate(
    Var1 = as.integer(
      sub(
        "non_norm_omega\\[(\\d+),\\d+\\]",
        "\\1",
        variable
      )
    ),
    Var2 = as.integer(
      sub(
        "non_norm_omega\\[\\d+,(\\d+)\\]",
        "\\1",
        variable
      )
    )
  ) %>%
  select(
    -variable
  )

ggplot(median_cor_matrix_df, aes(Var1, Var2)) +
  geom_point(
    aes(size = abs(Freq), fill = Freq),
    shape = 21,
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1),
    space = "Lab",
    name = "Correlation"
  ) +
  scale_size_continuous(
    range = c(1, 20),
    guide = "none"
  ) +
  coord_fixed() +
  ylab("") +
  xlab("") +
  ggtitle("Exponential Correlation Matrix Visual") +
  theme( # //axis.title.x = element_blank(),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 14),
    axis.line = element_line(colour = "black", size = 1.25),
    axis.ticks = element_line(colour = "black", size = 1.5),
    axis.ticks.length = unit(.25, "cm"),
    panel.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(colour = "white", face = "bold", size = 12),
    legend.text = element_text(colour = "white", face = "bold", size = 12),
    legend.background = element_rect(fill = "black"),
    legend.key.width = unit(.025, "npc"),
    plot.title = element_text(face = "bold", size = 16),
    strip.text = element_text(colour = "white", face = "bold", size = 12),
    strip.background = element_rect(fill = "black")
  )


cor_matrix_draws_df <- as.table(
  cor_matrix_draws
) %>%
  as.data.frame() %>%
  filter(
    chain == 1
  ) %>%
  select(
    -chain
  ) %>%
  mutate(
    Var1 = as.integer(
      sub(
        "non_norm_omega\\[(\\d+),\\d+\\]",
        "\\1",
        variable
      )
    ),
    Var2 = as.integer(
      sub(
        "non_norm_omega\\[\\d+,(\\d+)\\]",
        "\\1",
        variable
      )
    )
  ) %>%
  select(
    -variable
  ) %>%
  mutate(
    Var1 = case_when(
      Var1 == 1 ~ "Site 1",
      Var1 == 2 ~ "Site 2",
      Var1 == 3 ~ "Site 3",
      Var1 == 4 ~ "Site 4"
    ),
    Var2 = case_when(
      Var2 == 1 ~ "Site 1",
      Var2 == 2 ~ "Site 2",
      Var2 == 3 ~ "Site 3",
      Var2 == 4 ~ "Site 4"
    )
  ) %>%
  left_join(
    corr_df,
    by = c(
      "Var1",
      "Var2"
    )
  ) %>%
  mutate(
    Var1 = factor(
      Var1,
      levels = c(
        "Site 4",
        "Site 3",
        "Site 2",
        "Site 1"
      )
    ),
    Var2 = factor(
      Var2
    )
  )
ggplot(cor_matrix_draws_df) +
  geom_density(
    aes(x = Freq),
    fill = "blue",
    alpha = 0.5
  ) +
  geom_vline(
    aes(xintercept = value),
    color = "red",
    linetype = "dashed",
    size = 1.25
  ) +
  facet_grid(
    Var1 ~ Var2
  ) +
  ylab("Frequency") +
  xlab("Correlation") +
  ggtitle(
    "Stan Inf. Correlation Matrix W/out Dist. Visual (Red line is actual)"
  ) +
  theme( # //axis.title.x = element_blank(),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 11),
    axis.line = element_line(colour = "black", size = 1.25),
    axis.ticks = element_line(colour = "black", size = 1.5),
    axis.ticks.length = unit(.25, "cm"),
    panel.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(colour = "white", face = "bold", size = 12),
    legend.text = element_text(colour = "white", face = "bold", size = 12),
    legend.background = element_rect(fill = "black"),
    legend.key.width = unit(.025, "npc"),
    plot.title = element_text(face = "bold", size = 16),
    strip.text = element_text(colour = "white", face = "bold", size = 12),
    strip.background = element_rect(fill = "black")
  )
