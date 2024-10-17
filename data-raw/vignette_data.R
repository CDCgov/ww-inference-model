set.seed(1)
simulated_data <- wwinference::generate_simulated_data()
hosp_data_from_sim <- simulated_data$hosp_data
ww_data <- simulated_data$ww_data
ww_data_eval <- simulated_data$ww_data_eval
hosp_data <- hosp_data_from_sim |>
  dplyr::mutate("location" = "example state")
hosp_data_eval <- simulated_data$hosp_data_eval
subpop_hosp_data <- simulated_data$subpop_hosp_data
subpop_hosp_data_eval <- simulated_data$subpop_hosp_data_eval
true_global_rt <- simulated_data$true_global_rt

usethis::use_data(hosp_data, overwrite = TRUE)
usethis::use_data(hosp_data_eval, overwrite = TRUE)
usethis::use_data(ww_data, overwrite = TRUE)
usethis::use_data(ww_data_eval, overwrite = TRUE)
usethis::use_data(subpop_hosp_data, overwrite = TRUE)
usethis::use_data(subpop_hosp_data_eval, overwrite = TRUE)
usethis::use_data(true_global_rt, overwrite = TRUE)
