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
true_global_rt <- simulated_data$true_global_rt

usethis::use_data(hosp_data, overwrite = TRUE)
usethis::use_data(hosp_data_eval, overwrite = TRUE)
usethis::use_data(ww_data, overwrite = TRUE)
usethis::use_data(true_global_rt, overwrite = TRUE)
