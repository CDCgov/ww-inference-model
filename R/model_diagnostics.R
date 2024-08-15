#' Get a table of diagnostic flags
#'
#' @description
#' This function takes in the output from a cmdstanr$sample() function (the
#' fit object) and a series of diagnostic tolerances and returns
#' a dataframe containing flags for whether any of the diagnostic thresholds
#' were exceeded, which would indicate that the model did not properly
#' converge. This funtion has a default method that takes
#' the CmdStan fitting object, as well as an S3 method for objects of class
#' 'wwinference_fit'
#'
#'
#' @param x Either an object of the 'wwinference_fit' class or
#' the R6 Cmdstan Object fit object
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
#' @family diagnostics
#' @return flag_df: dataframe containing columns for each of the flags,
#' if any flags are TRUE that indicates some model issue
#' @export
get_model_diagnostic_flags <- function(x, ...) {
  UseMethod("get_model_diagnostic_flags")
}

#' S3 method for getting a table of diagnostic flags fpr a wwinference_fit
#' object
#'
#' This method overloads the generic get_model_diagnostic_flags function
#' specifically for objects of type 'wwinference_fit'.
#'
#' @rdname get_model_diagnostic_flags
#' @export
get_model_diagnostic_flags.wwinference_fit <- function(x, ...) {
  get_model_diagnostic_flags.default(
    x = x$fit$result
  )
}

#' @rdname get_model_diagnostic_flags
#' @export
get_model_diagnostic_flags.default <- function(x,
                                               ebmfi_tolerance = 0.2,
                                               divergences_tolerance = 0.01,
                                               frac_high_rhat_tolerance = 0.05,
                                               rhat_tolerance = 1.05,
                                               max_tree_depth_tol = 0.01) {
  n_chains <- x$num_chains()
  iter_sampling <- x$metadata()$iter_sampling

  # Summary is a large dataframe with diagnostics for each parameters
  summary <- x$summary()
  diagnostic_summary <- x$diagnostic_summary(quiet = TRUE)

  flag_low_embfi <- mean(diagnostic_summary$ebfmi) <= ebmfi_tolerance
  max_n_divergences <- n_chains * iter_sampling * divergences_tolerance
  flag_too_many_divergences <- any(
    diagnostic_summary$num_divergent >= max_n_divergences
  )
  frac_high_rhat <- as.numeric(mean(summary[, "rhat"]$rhat > rhat_tolerance,
    na.rm = TRUE
  ))
  flag_high_rhat <- frac_high_rhat >= frac_high_rhat_tolerance
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

  # Message if a flag doesn't pass. Still return
  # the same data, but we want user to know the issue
  if (any(flag_df[1, ])) {
    warning("Model flagged for convergence issues, run model diagnostics
      on the output stanfit object for more information")
  }
  return(flag_df)
}


#' Method for printing the CmdStan diagnostic summary
#' for a wwinference_fit object
#'
#' @param ww_fit An object of class wwinference_fit
#'
#' @family diagnostics
#' @export
diagnostic_summary <- function(ww_fit, ...) {
  ww_fit$fit$result$summary()
}

#' Method for printing the CmdStan parameter diagnostics for a
#' wwinference_fit_object
#'
#' @param ww_fit An object of class wwinference_fit
#'
#' @family diagnostics
#' @export
parameter_diagnostics <- function(ww_fit, ...) {
  ww_fit$fit$result$diagnostic_summary(quiet = TRUE)
}
