#### Final plots for spatial ww paper
### Run spatial_wwinference.Rmd then spatial_ww_paper_script.R(from scratch)
### before this script

## Simulation Presets w/ ww total plot
sim_presets_ww_plot <- (fake_locs_plot + central_rt_curve) /
  corr_matrix_plot /
  (ww_total_plot +
    labs(
      title = "Simulated Wastewater Data"
    )) +
  plot_layout(
    heights = c(1, 1.25, 1.25)
  ) &
  plot_annotation(tag_levels = "A")
sim_presets_ww_plot

## Simulated Data
sim_data_plot

## Fits Plot
hosp_ww_result
rt_results

## Hosp eval
# averaged across date
hosp_eval_total <- hosp_rel_eval_total_crps /
  (plot_spacer() + hosp_eval_total_alt + plot_spacer()) +
  plot_layout(
    heights = c(1, 4)
  ) &
  plot_annotation(
    tag_levels = "A"
  )
hosp_eval_total
# by date
hosp_eval_by_date <- hosp_eval_by_date_crps /
  hosp_eval_by_date_rel_crps /
  hosp_eval_by_date_bias +
  plot_layout(
    axis_titles = "collect"
  ) &
  plot_annotation(
    tag_levels = "A"
  )
hosp_eval_by_date
# quantile-quantile
hosp_coverage_plot

## WW Eval
# by site
ww_eval_by_site <- ww_eval_by_site_crps /
  ww_eval_by_site_rel_crps /
  ww_eval_by_site_bias +
  plot_annotation(
    tag_levels = "A"
  )
ww_eval_by_site
# by date
ww_eval_by_date <- ww_eval_by_date_crps /
  ww_eval_by_date_rel_crps /
  ww_eval_by_date_bias +
  plot_annotation(
    tag_levels = "A"
  )
ww_eval_by_date
# by date and site
ww_eval_by_site_date <- ww_eval_by_site_date_crps /
  ww_eval_by_site_date_rel_crps /
  ww_eval_by_site_date_bias +
  plot_annotation(
    tag_levels = "A"
  )
ww_eval_by_site_date
# quantile-quantile
ww_coverage_by_site
ww_coverage
