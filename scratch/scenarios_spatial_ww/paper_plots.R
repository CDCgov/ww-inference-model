# Manuscript plots
# RUN AFTER RUNNING SPATIAL VIGNETTE!!!
(
  fake_locs_plot +
    central_rt_curve +
    ggtitle("Central Reproduction Number") +
    ylab("R(t) Value")
) /
  corr_matrix_plot +
  theme(legend.position = "bottom")


(
  hosp_data_plot +
    ggtitle("Hospital Admissions") +
    global_rt_data_plot +
    ggtitle("Global Reproduction Number")
) / (
  ww_data_plot +
    ggtitle("Pathogen Genome Conc. in Wastewater per Site") +
    site_rt_data_plot +
    ggtitle("Site Reproduction Number")
)


ww_total_plot +
  ggtitle("Pathogen Genome Conc. in Wastewater per Site")


(
  hosp_result_plot +
    ggtitle("Hospital Admissions Inference")
) /
  (
    ww_result_plot +
      ggtitle("Pathogen Genome Conc. in Wastewater per Site Inference")
  ) +
  plot_layout(
    guides = "collect",
    heights = c(1, 2)
  ) &
  theme(legend.position = "bottom")


(
  global_rt_result_plot +
    ggtitle("Global Reproduction Number")
) /
  (
    site_rt_result_plot +
      ggtitle("Site Reproduction Number")
  ) +
  plot_layout(
    guides = "collect",
    heights = c(1, n_sites)
  ) &
  theme(legend.position = "bottom")


iid_corr_matrix_result_plot + exp_corr_matrix_result_plot +
  rand_corr_matrix_result_plot +
  plot_layout(
    guides = "collect"
  ) &
  theme(legend.position = "bottom")



(hosp_eval_nowcast_plot + hosp_eval_forecast_plot) /
  (hosp_by_date_eval_nowcast_plot + hosp_by_date_eval_forecast_plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
