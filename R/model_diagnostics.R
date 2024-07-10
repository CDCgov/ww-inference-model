#' Get diagnostic flags
#'
#' @description
#' This function takes in the output from a cmdstanr$sample() function (the
#' fit object) and a series of diagnostic tolerances and returns
#' a dataframe containing flags for whether any of the diagnostic thresholds
#' were exceeded, which would indicate that the model did not properly
#' converge
#'
#'
#' @param stan_fit_object The R6 Cmdstan Object fit object
#' @param ebmfi_tolerance float indicating the tolerance for EBMFI
#' (estimated bayesian fraction of missing information), default is `0.2`
#' @param divergences_tolerance float indicating the tolerance for the
#' proportion of sampling iterations that are divergent, default is `0.01`
#' @param frac_high_rhat_tolerance float indicating the tolerance for the
#' proportion of parameters rhats>rhat_tolderance, default is `0.05`
#' @param rhat_tolerance float indicating the tolerance for the rhat for
#' individual parameters, default is `1.05`
#' @param max_tree_depth_tol float indicating the tolerance for the proportion
#' of iterations that exceed the maximum tree depth, default is `0.01`
#'
#' @return flag_df: dataframe containing columns for each of the flags,
#' if any flags are TRUE that indicates some model issue
#' @export
#'
get_model_diagnostic_flags <- function(stan_fit_object,
                                       ebmfi_tolerance = 0.2,
                                       divergences_tolerance = 0.01,
                                       frac_high_rhat_tolerance = 0.05,
                                       rhat_tolerance = 1.05,
                                       max_tree_depth_tol = 0.01) {
  n_chains <- stan_fit_object$num_chains()
  iter_sampling <- stan_fit_object$metadata()$iter_sampling

  # Summary is a large dataframe with diagnostics for each parameters
  summary <- stan_fit_object$summary()
  diagnostic_summary <- stan_fit_object$diagnostic_summary(quiet = TRUE)

  flag_low_embfi <- mean(diagnostic_summary$ebfmi) <= ebmfi_tolerance
  max_n_divergences <- n_chains * iter_sampling * divergences_tolerance
  flag_too_many_divergences <- any(
    diagnostic_summary$num_divergent >= max_n_divergences
  )
  frac_high_rhat <- as.numeric(mean(summary[, "rhat"]$rhat > rhat_tolerance,
    na.rm = TRUE
  ))
  flag_high_rhat <- frac_high_rhat >=
    frac_high_rhat_tolerance
  max_n_max_treedepth <- n_chains * iter_sampling * max_tree_depth_tol
  flag_high_max_treedepth <- any(
    diagnostic_summary$num_max_tree_depth >= max_n_max_treedepth
  )

  flag_df <- data.frame(
    flag_high_max_treedepth,
    flag_too_many_divergences,
    flag_high_rhat,
    flag_low_embfi
  )
  return(flag_df)
}
