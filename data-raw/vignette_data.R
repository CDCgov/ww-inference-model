set.seed(1)
simulated_data <- wwinference::generate_simulated_data()
hosp_data <- simulated_data$hosp_data
ww_data <- simulated_data$ww_data
hosp_data_eval <- simulated_data$hosp_data_eval

usethis::use_data(hosp_data, overwrite = TRUE)
usethis::use_data(hosp_data_eval, overwrite = TRUE)
usethis::use_data(ww_data, overwrite = TRUE)
