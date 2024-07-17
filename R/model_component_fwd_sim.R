get_global_rt <- function(r_in_weeks,
                          n_weeks,
                          n_days,
                          global_rt_sd) {
  # Generate the state level weekly R(t) before infection feedback
  # Adds a bit of noise, can add more...
  r_weeks_w_noise <- (r_in_weeks * rnorm(
    length(r_in_weeks),
    1, global_rt_sd
  ))[1:n_weeks]

  # Convert to daily for input into renewal equation
  ind_m <- get_ind_m(n_days, n_weeks)
  r_days <- ind_m %*% r_weeks_w_noise
  return(r_days)
}
