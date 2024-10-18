set.seed(1)
simulated_data <- wwinference::generate_simulated_data()
hosp_data_from_sim <- simulated_data$hosp_data
ww_data_from_sim <- simulated_data$ww_data
# Add some columns and reorder sites to ensure package works as expected
# even if sites are not in order
ww_data <- ww_data_from_sim |>
  dplyr::mutate(
    "location" = "example state",
    "site" = .data$site + 1
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(desc(.data$site))
hosp_data <- hosp_data_from_sim |>
  dplyr::mutate("location" = "example state")
hosp_data_eval <- simulated_data$hosp_data_eval
rt_site_data <- simulated_data$rt_site_data
rt_global_data <- simulated_data$rt_global_data


set.seed(1)
simulated_data_ind <- wwinference::generate_simulated_data(
  use_spatial_corr = FALSE,
  aux_site_bool = FALSE
)
hosp_data_ind <- simulated_data_ind$hosp_data
ww_data_ind <- simulated_data_ind$ww_data
hosp_data_eval_ind <- simulated_data_ind$hosp_data_eval
rt_site_data_ind <- simulated_data_ind$rt_site_data
rt_global_data_ind <- simulated_data_ind$rt_global_data

usethis::use_data(hosp_data, overwrite = TRUE)
usethis::use_data(hosp_data_eval, overwrite = TRUE)
usethis::use_data(ww_data, overwrite = TRUE)
usethis::use_data(rt_site_data, overwrite = TRUE)
usethis::use_data(rt_global_data, overwrite = TRUE)


usethis::use_data(hosp_data_ind, overwrite = TRUE)
usethis::use_data(hosp_data_eval_ind, overwrite = TRUE)
usethis::use_data(ww_data_ind, overwrite = TRUE)
usethis::use_data(rt_site_data_ind, overwrite = TRUE)
usethis::use_data(rt_global_data_ind, overwrite = TRUE)
