## Run after spatial_wwinfernce markdown is ran
## Makes plots for spatial ww paper


# Simulation Presets
sim_presets_plot <- (fake_locs_plot + central_rt_curve) /
  (corr_matrix_plot) +
  plot_layout(
    heights = c(1, 1.5)
  ) &
  plot_annotation(
    tag_levels = "A"
  )
sim_presets_plot

# Simulated Data
rt_data_plot <- global_rt_data_plot / site_rt_data_plot +
  plot_layout(
    heights = c(1, 5)
  )
sim_data_plot <- ((hosp_data_plot / ww_data_plot) | rt_data_plot) +
  plot_annotation(
    tag_levels = "A"
  )
sim_data_plot

# Simulated Data By Color
ww_total_plot

# Hosp WW Inference
hosp_ww_result <- hosp_result_plot / ww_result_plot +
  plot_annotation(
    tag_levels = "A"
  ) +
  plot_layout(
    guides = "collect",
    heights = c(1, 2)
  ) &
  theme(
    legend.position = "bottom"
  )
hosp_ww_result

# Rt Inference
rt_results <- global_rt_result_plot / site_rt_result_plot +
  plot_layout(
    guides = "collect",
    heights = c(1, 5)
  ) +
  plot_annotation(
    tag_levels = "A"
  )
rt_results

# Correlations Results
corr_results <- iid_corr_matrix_result_plot + exp_corr_matrix_result_plot +
  rand_corr_matrix_result_plot +
  plot_annotation(
    tag_levels = "A"
  ) +
  plot_layout(
    guides = "collect"
  ) &
  theme(legend.position = "bottom")
corr_results

# Hosp Eval Scores
hosp_obj_for_eval <- all_pred_draws_df %>%
  filter(
    date > max(hosp_data$date)
  ) %>%
  inner_join(
    hosp_data_eval_full %>%
      filter(
        date > max(hosp_data$date)
      ) %>%
      mutate(
        gen_model_type = model_type,
        true_value = daily_hosp_admits_for_eval
      ) %>%
      select(
        date,
        true_value,
        gen_model_type
      ),
    by = c("date", "gen_model_type")
  ) %>%
  mutate(
    predicted = pred_value,
    observed = true_value,
    sample_id = draw,
    model = inf_model_type
  ) %>%
  select(
    predicted,
    observed,
    sample_id,
    date,
    model,
    gen_model_type
  )
hosp_scores <- hosp_obj_for_eval %>%
  as_forecast_sample(
    forecast_unit = c(
      "date", "model", "gen_model_type"
    )
  ) %>%
  score() %>%
  select(
    date, model, gen_model_type, bias, crps
  )
hosp_scores_overall <- hosp_scores %>%
  group_by(
    model, gen_model_type
  ) %>%
  summarise(
    bias = mean(bias),
    crps = mean(crps)
  ) %>%
  ungroup()
hosp_scores_overall_rel <- rbind(
  hosp_scores_overall %>%
    filter(
      gen_model_type == "IID"
    ) %>%
    mutate(
      Bias = bias / bias[which(model == "IID")],
      CRPS = crps / crps[which(model == "IID")]
    ),
  hosp_scores_overall %>%
    filter(
      gen_model_type == "Exponential"
    ) %>%
    mutate(
      Bias = bias / bias[which(model == "IID")],
      CRPS = crps / crps[which(model == "IID")]
    ),
  hosp_scores_overall %>%
    filter(
      gen_model_type == "Pre-set Corr. Matrix"
    ) %>%
    mutate(
      Bias = bias / bias[which(model == "IID")],
      CRPS = crps / crps[which(model == "IID")]
    )
) %>%
  filter(
    model != "IID"
  ) %>%
  select(
    -c(bias, crps)
  ) %>%
  mutate(
    Bias = abs(Bias)
  ) %>%
  pivot_longer(
    -c(model, gen_model_type),
    names_to = "Metric",
    values_to = "Value"
  )

hosp_rel_eval_total_crps <- ggplot(
  hosp_scores_overall_rel %>%
    filter(Metric != "Bias")
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  geom_point(
    aes(
      x = Value,
      y = Metric,
      shape = model,
      colour = model
    ),
    size = 2.5
  ) +
  scale_colour_manual(
    breaks = c("Exponential", "Unstructured", "IID"),
    values = c("darkviolet", "deeppink3", "darksalmon")
  ) +
  scale_shape_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c(6, 0, 8)
  ) +
  facet_grid(~gen_model_type) +
  scale_x_continuous(
    trans = "log10",
    breaks = c(0.75, 0.9, 1, 1.1, 1.25)
  ) +
  coord_cartesian(
    xlim = c(0.75, 1.325)
  ) +
  labs(
    shape = "Assumed Correlation Structure",
    color = "Assumed Correlation Structure",
    x = "Relative Score"
  ) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    legend.position = "bottom"
  )
hosp_rel_eval_total_crps


hosp_eval_total_bias <- ggplot(
  hosp_scores_overall %>%
    pivot_longer(
      -c(model, gen_model_type),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(
      Metric = case_when(
        Metric == "bias" ~ "Bias",
        .default = Metric
      )
    ) %>%
    filter(
      Metric == "Bias"
    )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    aes(
      x = Value,
      y = Metric,
      shape = model
    )
  ) +
  scale_shape_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c(6, 0, 8)
  ) +
  facet_grid(~gen_model_type) +
  labs(
    shape = "Assumed Correlation Structure",
    x = "Bias"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title.y = element_blank()
  )
hosp_eval_total_bias
hosp_eval_total <- (hosp_rel_eval_total_crps / hosp_eval_total_bias) +
  plot_layout(
    axis_titles = "collect"
  ) &
  plot_annotation(
    tag_levels = "A"
  )
hosp_eval_total


hosp_eval_total_alt <- ggplot(
  hosp_scores_overall %>%
    pivot_longer(
      -c(model, gen_model_type),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(
      Metric = case_when(
        Metric == "bias" ~ "Bias",
        Metric == "crps" ~ "CRPS",
        .default = Metric
      )
    )
) +
  geom_col(
    aes(
      x = gen_model_type,
      y = Value,
      colour = model,
      fill = model
    ),
    alpha = 0.25,
    position = "dodge"
  ) +
  facet_grid2(
    Metric ~ .,
    scales = "free"
  ) +
  scale_colour_manual(
    breaks = c("Exponential", "Unstructured", "IID"),
    values = c("darkviolet", "deeppink3", "darksalmon")
  ) +
  scale_fill_manual(
    breaks = c("Exponential", "Unstructured", "IID"),
    values = c("darkviolet", "deeppink3", "darksalmon")
  ) +
  labs(
    color = "Assumed Correlation Structure",
    fill = "Assumed Correlation Structure",
    y = "Score",
    x = "Actual Corr. Structure"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )
hosp_eval_total_alt


hosp_score_by_date_rel <- rbind(
  hosp_scores %>%
    filter(
      gen_model_type == "IID"
    ) %>%
    group_by(
      date
    ) %>%
    mutate(
      `Relative CRPS` = crps / crps[which(model == "IID")]
    ),
  hosp_scores %>%
    filter(
      gen_model_type == "Exponential"
    ) %>%
    group_by(
      date
    ) %>%
    mutate(
      `Relative CRPS` = crps / crps[which(model == "IID")]
    ),
  hosp_scores %>%
    filter(
      gen_model_type == "Pre-set Corr. Matrix"
    ) %>%
    group_by(
      date
    ) %>%
    mutate(
      `Relative CRPS` = crps / crps[which(model == "IID")]
    )
) %>%
  mutate(
    Bias = bias,
    CRPS = crps
  ) %>%
  select(
    -c(bias, crps)
  ) %>%
  pivot_longer(
    -c(date, model, gen_model_type),
    names_to = "Metric",
    values_to = "Value"
  )
hosp_eval_by_date_bias <- ggplot(
  hosp_score_by_date_rel %>%
    filter(
      Metric == "Bias"
    )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = as.Date(forecast_date)
  ) +
  geom_line(
    aes(
      x = date,
      y = Value,
      colour = model
    )
  ) +
  geom_point(
    aes(
      x = date,
      y = Value,
      colour = model,
      shape = model
    ),
    size = 2.25
  ) +
  scale_shape_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c(6, 0, 8)
  ) +
  scale_colour_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c("darksalmon", "darkviolet", "deeppink3")
  ) +
  facet_grid2(
    Metric ~ gen_model_type,
    scales = "free",
    independent = "y"
  ) +
  labs(
    color = "Assumed Corr. Structure",
    shape = "Assumed Corr. Structure",
    y = "Bias",
    x = "Date"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )
hosp_eval_by_date_rel_crps <- ggplot(
  hosp_score_by_date_rel %>%
    filter(
      Metric == "Relative CRPS",
      model != "IID"
    )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = as.Date(forecast_date)
  ) +
  geom_line(
    aes(
      x = date,
      y = Value,
      colour = model
    )
  ) +
  geom_point(
    aes(
      x = date,
      y = Value,
      colour = model,
      shape = model
    ),
    size = 2.25
  ) +
  scale_y_continuous(
    trans = "log10",
    breaks = c(0.5, 0.7, 1, 1.5, 2)
  ) +
  coord_cartesian(
    ylim = c(0.4, 2.5)
  ) +
  scale_shape_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c(6, 0, 8)
  ) +
  scale_colour_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c("darksalmon", "darkviolet", "deeppink3")
  ) +
  facet_grid2(
    Metric ~ gen_model_type
  ) +
  labs(
    color = "Assumed Corr. Structure",
    shape = "Assumed Corr. Structure",
    y = "Relative CRPS",
    x = "Date"
  ) +
  theme_bw() +
  guides(
    color = "none",
    shape = "none"
  )
hosp_eval_by_date_rel_crps
hosp_eval_by_date_crps <- ggplot(
  hosp_score_by_date_rel %>%
    filter(
      Metric == "CRPS"
    )
) +
  geom_vline(
    xintercept = as.Date(forecast_date)
  ) +
  geom_line(
    aes(
      x = date,
      y = Value,
      colour = model
    )
  ) +
  geom_point(
    aes(
      x = date,
      y = Value,
      colour = model,
      shape = model
    ),
    size = 2.25
  ) +
  scale_shape_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c(6, 0, 8)
  ) +
  scale_colour_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c("darksalmon", "darkviolet", "deeppink3")
  ) +
  facet_grid2(
    Metric ~ gen_model_type,
    scales = "free",
    independent = "y"
  ) +
  labs(
    color = "Assumed Corr. Structure",
    shape = "Assumed Corr. Structure",
    y = "Actual CRPS",
    x = "Date"
  ) +
  theme_bw() +
  guides(
    color = "none",
    shape = "none"
  )
hosp_eval_by_date_crps


hosp_obj_for_eval <- hosp_obj_for_eval %>%
  as_forecast_sample(
    forecast_unit = c(
      "date", "model", "gen_model_type"
    )
  )
hosp_obj_for_quant_eval <- as_forecast_quantile(
  hosp_obj_for_eval,
  probs = c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99),
  type = 7,
  forecast_unit = c(
    "date", "model", "gen_model_type"
  )
)
hosp_coverage_obj <- hosp_obj_for_quant_eval %>%
  get_coverage(
    by = c(
      "model", "gen_model_type"
    )
  )
hosp_coverage_plot <- ggplot(
  hosp_coverage_obj
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  geom_line(
    aes(
      x = quantile_level,
      y = quantile_coverage,
      colour = model
    )
  ) +
  scale_colour_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c("darksalmon", "darkviolet", "deeppink3")
  ) +
  facet_grid2(~gen_model_type) +
  labs(
    color = "Assumed Correlation Structure",
    y = "% Obs below quantile level",
    x = "Quantile level"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )




# WW Eval Scores
ww_obj_for_eval_forecast <- all_ww_draws_df %>%
  filter(
    date > max(hosp_data$date)
  ) %>%
  inner_join(
    ww_data_eval_full %>%
      filter(
        date > max(hosp_data$date)
      ) %>%
      mutate(
        gen_model_type = model_type,
        true_value = log_genome_copies_per_ml_eval,
        subpop_name = case_when(
          site == 1 ~ "Site: 1",
          site == 2 ~ "Site: 2",
          site == 3 ~ "Site: 3",
          site == 4 ~ "Site: 4"
        )
      ) %>%
      select(
        date,
        true_value,
        subpop_name,
        gen_model_type
      ),
    by = c("date", "subpop_name", "gen_model_type")
  ) %>%
  mutate(
    predicted = pred_value,
    observed = true_value,
    sample_id = draw,
    model = inf_model_type
  ) %>%
  select(
    predicted,
    observed,
    sample_id,
    date,
    subpop_name,
    model,
    gen_model_type
  )
ww_scores <- ww_obj_for_eval_forecast %>%
  as_forecast_sample(
    forecast_unit = c(
      "date", "subpop_name", "model", "gen_model_type"
    )
  ) %>%
  score() %>%
  select(
    date, subpop_name, model, gen_model_type, bias, crps
  )
ww_score_by_date_site <- rbind(
  ww_scores %>%
    filter(
      gen_model_type == "IID"
    ) %>%
    group_by(
      date, subpop_name
    ) %>%
    mutate(
      `Relative CRPS` = crps / crps[which(model == "IID")]
    ),
  ww_scores %>%
    filter(
      gen_model_type == "Exponential"
    ) %>%
    group_by(
      date, subpop_name
    ) %>%
    mutate(
      `Relative CRPS` = crps / crps[which(model == "IID")]
    ),
  ww_scores %>%
    filter(
      gen_model_type == "Pre-set Corr. Matrix"
    ) %>%
    group_by(
      date, subpop_name
    ) %>%
    mutate(
      `Relative CRPS` = crps / crps[which(model == "IID")]
    )
) %>%
  mutate(
    Bias = bias,
    CRPS = crps
  ) %>%
  select(
    -c(bias, crps)
  ) %>%
  pivot_longer(
    -c(date, model, gen_model_type, subpop_name),
    names_to = "Metric",
    values_to = "Value"
  )
ww_score_by_date_site_wovall <- rbind(
  ww_score_by_date_site,
  ww_score_by_date_site %>%
    group_by(date, model, gen_model_type, Metric) %>%
    summarize(
      subpop_name = "Overall",
      Value = mean(Value, na.rm = TRUE),
      .groups = "drop"
    )
) %>%
  mutate(
    subpop_name = factor(subpop_name)
  )

ww_eval_by_site_rel_crps <- ggplot(
  ww_score_by_date_site_wovall %>%
    group_by(
      subpop_name, model, gen_model_type, Metric
    ) %>%
    summarise(
      Value = mean(Value),
      .groups = "drop"
    ) %>%
    filter(
      Metric == "Relative CRPS",
      model != "IID"
    ) %>%
    mutate(
      subpop_name = factor(
        subpop_name,
        levels = c("Site: 4", "Site: 3", "Site: 2", "Site: 1", "Overall")
      )
    )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  geom_point(
    aes(
      x = Value,
      y = subpop_name,
      shape = model,
      colour = model
    ),
    size = 2.25
  ) +
  scale_x_continuous(
    trans = "log10",
    breaks = c(0.5, 0.7, 1, 1.5, 2)
  ) +
  coord_cartesian(
    xlim = c(0.4, 2.5)
  ) +
  scale_colour_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c("darksalmon", "darkviolet", "deeppink3")
  ) +
  scale_shape_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c(6, 0, 8)
  ) +
  facet_grid2(
    Metric ~ gen_model_type,
    scales = "free",
    independent = "x"
  ) +
  labs(
    shape = "Assumed Correlation Structure",
    color = "Assumed Correlation Structure",
    x = "Relative CRPS"
  ) +
  guides(
    shape = "none",
    color = "none"
  ) +
  theme_bw() +
  theme(
    axis.title.y = element_blank()
  )
ww_eval_by_site_crps <- ggplot(
  ww_score_by_date_site_wovall %>%
    group_by(
      subpop_name, model, gen_model_type, Metric
    ) %>%
    summarise(
      Value = mean(Value),
      .groups = "drop"
    ) %>%
    filter(
      Metric == "CRPS"
    ) %>%
    mutate(
      subpop_name = factor(
        subpop_name,
        levels = c("Site: 4", "Site: 3", "Site: 2", "Site: 1", "Overall")
      )
    )
) +
  geom_point(
    aes(
      x = Value,
      y = subpop_name,
      shape = model,
      colour = model
    ),
    size = 2.25
  ) +
  scale_colour_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c("darksalmon", "darkviolet", "deeppink3")
  ) +
  scale_shape_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c(6, 0, 8)
  ) +
  facet_grid2(
    Metric ~ gen_model_type
  ) +
  labs(
    shape = "Assumed Correlation Structure",
    color = "Assumed Correlation Structure",
    x = "CRPS"
  ) +
  guides(
    shape = "none",
    color = "none"
  ) +
  theme_bw() +
  theme(
    axis.title.y = element_blank()
  )
ww_eval_by_site_bias <- ggplot(
  ww_score_by_date_site_wovall %>%
    group_by(
      subpop_name, model, gen_model_type, Metric
    ) %>%
    summarise(
      Value = mean(Value),
      .groups = "drop"
    ) %>%
    filter(
      Metric == "Bias"
    ) %>%
    mutate(
      subpop_name = factor(
        subpop_name,
        levels = c("Site: 4", "Site: 3", "Site: 2", "Site: 1", "Overall")
      )
    )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    aes(
      x = Value,
      y = subpop_name,
      shape = model,
      colour = model
    ),
    size = 2.25
  ) +
  scale_colour_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c("darksalmon", "darkviolet", "deeppink3")
  ) +
  scale_shape_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c(6, 0, 8)
  ) +
  facet_grid2(
    Metric ~ gen_model_type
  ) +
  labs(
    shape = "Assumed Correlation Structure",
    color = "Assumed Correlation Structure",
    x = "Bias"
  ) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    legend.position = "bottom"
  )



ww_eval_by_date_rel_crps <- ggplot(
  ww_score_by_date_site %>%
    group_by(
      date, model, gen_model_type, Metric
    ) %>%
    summarise(
      Value = mean(Value)
    ) %>%
    filter(
      Metric == "Relative CRPS",
      model != "IID"
    )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = as.Date(forecast_date)
  ) +
  geom_line(
    aes(
      x = date,
      y = Value,
      colour = model
    )
  ) +
  geom_point(
    aes(
      x = date,
      y = Value,
      colour = model,
      shape = model
    ),
    size = 2.25
  ) +
  scale_y_continuous(
    trans = "log10",
    breaks = c(0.5, 0.7, 1, 1.5, 2)
  ) +
  coord_cartesian(
    ylim = c(0.5, 2)
  ) +
  scale_shape_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c(6, 0, 8)
  ) +
  scale_colour_manual(
    breaks = c("Exponential", "Unstructured", "IID"),
    values = c("darkviolet", "deeppink3", "darksalmon")
  ) +
  facet_grid2(
    Metric ~ gen_model_type,
    scales = "free",
    independent = "y"
  ) +
  labs(
    color = "Assumed Correlation Structure",
    shape = "Assumed Correlation Structure",
    y = "Relative CRPS",
    x = "Date"
  ) +
  guides(
    shape = "none",
    color = "none"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )
ww_eval_by_date_crps <- ggplot(
  ww_score_by_date_site %>%
    group_by(
      date, model, gen_model_type, Metric
    ) %>%
    summarise(
      Value = mean(Value)
    ) %>%
    filter(
      Metric == "CRPS"
    )
) +
  geom_vline(
    xintercept = as.Date(forecast_date)
  ) +
  geom_line(
    aes(
      x = date,
      y = Value,
      colour = model
    )
  ) +
  geom_point(
    aes(
      x = date,
      y = Value,
      colour = model,
      shape = model
    ),
    size = 2.25
  ) +
  scale_shape_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c(6, 0, 8)
  ) +
  scale_colour_manual(
    breaks = c("Exponential", "Unstructured", "IID"),
    values = c("darkviolet", "deeppink3", "darksalmon")
  ) +
  facet_grid2(
    Metric ~ gen_model_type
  ) +
  labs(
    color = "Assumed Correlation Structure",
    shape = "Assumed Correlation Structure",
    y = "CRPS",
    x = "Date"
  ) +
  guides(
    shape = "none",
    color = "none"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )
ww_eval_by_date_bias <- ggplot(
  ww_score_by_date_site %>%
    group_by(
      date, model, gen_model_type, Metric
    ) %>%
    summarise(
      Value = mean(Value)
    ) %>%
    filter(
      Metric == "Bias"
    )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = as.Date(forecast_date)
  ) +
  geom_line(
    aes(
      x = date,
      y = Value,
      colour = model
    )
  ) +
  geom_point(
    aes(
      x = date,
      y = Value,
      colour = model,
      shape = model
    ),
    size = 2.25
  ) +
  scale_shape_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c(6, 0, 8)
  ) +
  scale_colour_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c("darksalmon", "darkviolet", "deeppink3")
  ) +
  facet_grid2(
    Metric ~ gen_model_type
  ) +
  labs(
    color = "Assumed Correlation Structure",
    shape = "Assumed Correlation Structure",
    y = "Bias",
    x = "Date"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )



ww_eval_by_site_date_rel_crps <- ggplot(
  ww_score_by_date_site_wovall %>%
    filter(
      Metric == "Relative CRPS",
      model != "IID"
    )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = as.Date(forecast_date)
  ) +
  geom_line(
    aes(
      x = date,
      y = Value,
      colour = subpop_name,
      linetype = model
    )
  ) +
  geom_point(
    aes(
      x = date,
      y = Value,
      colour = subpop_name,
      shape = model
    ),
    size = 2.25
  ) +
  facet_grid2(
    Metric ~ gen_model_type,
    scales = "free",
    independent = "y"
  ) +
  scale_y_continuous(
    trans = "log10",
    breaks = c(0.5, 0.7, 1, 1.5, 2)
  ) +
  coord_cartesian(
    ylim = c(1 / 2.75, 2.75)
  ) +
  scale_shape_manual(
    breaks = c("Exponential", "Unstructured", "IID"),
    values = c(0, 8, 6)
  ) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    colour = "Site",
    linetype = "Assumed Correlation Structure",
    shape = "Assumed Correlation Structure",
    y = "Relative CRPS",
    x = "Date"
  ) +
  guides(
    color = "none",
    linetype = "none",
    shape = "none"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )
ww_eval_by_site_date_crps <- ggplot(
  ww_score_by_date_site_wovall %>%
    filter(
      Metric == "CRPS"
    )
) +
  geom_vline(
    xintercept = as.Date(forecast_date)
  ) +
  geom_line(
    aes(
      x = date,
      y = Value,
      colour = subpop_name,
      linetype = model
    )
  ) +
  geom_point(
    aes(
      x = date,
      y = Value,
      colour = subpop_name,
      shape = model
    ),
    size = 2.25
  ) +
  facet_grid2(
    Metric ~ gen_model_type
  ) +
  scale_shape_manual(
    breaks = c("Exponential", "Unstructured", "IID"),
    values = c(0, 8, 6)
  ) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    colour = "Site",
    linetype = "Assumed Correlation Structure",
    shape = "Assumed Correlation Structure",
    y = "CRPS",
    x = "Date"
  ) +
  guides(
    color = "none",
    linetype = "none",
    shape = "none"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )
ww_eval_by_site_date_bias <- ggplot(
  ww_score_by_date_site_wovall %>%
    filter(
      Metric == "Bias"
    )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = as.Date(forecast_date)
  ) +
  geom_line(
    aes(
      x = date,
      y = Value,
      colour = subpop_name,
      linetype = model
    )
  ) +
  geom_point(
    aes(
      x = date,
      y = Value,
      colour = subpop_name,
      shape = model
    ),
    size = 2.25
  ) +
  facet_grid2(
    Metric ~ gen_model_type
  ) +
  scale_shape_manual(
    breaks = c("Exponential", "Unstructured", "IID"),
    values = c(0, 8, 6)
  ) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    colour = "Site",
    linetype = "Assumed Correlation Structure",
    shape = "Assumed Correlation Structure",
    y = "Bias",
    x = "Date"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )




ww_obj_for_eval_forecast <- ww_obj_for_eval_forecast %>%
  as_forecast_sample(
    forecast_unit = c(
      "date", "subpop_name", "model", "gen_model_type"
    )
  )
ww_obj_for_quant_eval <- as_forecast_quantile(
  ww_obj_for_eval_forecast,
  probs = c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99),
  type = 7,
  forecast_unit = c(
    "date", "subpop_name", "model", "gen_model_type"
  )
)
ww_coverage_obj_by_site <- ww_obj_for_quant_eval %>%
  get_coverage(
    by = c(
      "subpop_name", "model", "gen_model_type"
    )
  )
ww_coverage_obj <- ww_obj_for_quant_eval %>%
  get_coverage(
    by = c(
      "model", "gen_model_type"
    )
  )
ww_coverage_by_site <- ggplot(
  ww_coverage_obj_by_site
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  geom_line(
    aes(
      x = quantile_level,
      y = quantile_coverage,
      colour = model
    )
  ) +
  scale_colour_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c("darksalmon", "darkviolet", "deeppink3")
  ) +
  facet_grid2(subpop_name ~ gen_model_type) +
  labs(
    color = "Assumed Correlation Structure",
    y = "% Obs below quantile level",
    x = "Quantile level"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

ww_coverage <- ggplot(
  ww_coverage_obj
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  geom_line(
    aes(
      x = quantile_level,
      y = quantile_coverage,
      colour = model
    )
  ) +
  scale_colour_manual(
    breaks = c("IID", "Exponential", "Unstructured"),
    values = c("darksalmon", "darkviolet", "deeppink3")
  ) +
  facet_grid2(~gen_model_type) +
  labs(
    color = "Assumed Correlation Structure",
    y = "% Obs below quantile level",
    x = "Quantile level"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )
