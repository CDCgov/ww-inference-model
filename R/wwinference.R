wwinference <- function(ww_data,
                        count_data,
                        model_spec = get_model_spec(
                          forecast_date =
                            "2023-12-06"
                        ),
                        mcmc_options = get_mcmc_options(),
                        compiled_model = compile_model()) {
  # Check that data is compatible with specifications
  check_date(ww_data, model_spec$forecast_date)
  check_date(count_data, model_spec$forecast_date)

  # If checks pass, create stan data object
  stan_data <- get_stan_data(
    input_count_data = count_data,
    input_ww_data = ww_data,
    forecast_date = model_spec$forecast_date,
    calibration_time = model_spec$calibration_time,
    forecast_horizon = model_spec$forecast_horizon,
    generation_interval = model_spec$generation_interval,
    inf_to_count_delay = model_spec$inf_to_count_delay,
    infection_feedback_pmf = model_spec$infection_feedback_pmf,
    params = model_spec$params,
    compute_likelihood = 1
  )

  init_lists <- c()
  for (i in 1:model_spec$n_chains) {
    init_lists[[i]] <- get_inits(stan_data, params)
  }


  fit_model <- function(compiled_model,
                        standata,
                        model_spec,
                        init_lists) {
    fit <- compiled_model$sample(
      data = stan_data,
      init = init_lists,
      seed = model_spec$seed,
      iter_sampling = model_spec$iter_sampling,
      iter_warmup = model_spec$iter_warmup,
      max_treedepth = model_spec$max_treedepth,
      chains = model_spec$n_chains,
      parallel_chains = model_spec$n_chains
    )
    print(fit)
    return(fit)
  }

  # This returns the cmdstan object if the model runs, and result = NULL if
  # the model errors
  safe_fit_model <- purrr::safely(fit_model)

  fit <- safe_fit_model(
    compiled_model,
    standata,
    model_spec,
    init_lists
  )

  if (!is.null(fit$error)) { # If the model errors, return a list with the
    # error and everything else NULL
    out <- list(
      error = fit$error[[1]]
    )
    message(error)
  } else {
    draws <- fit$result$draws()
    diagnostics <- fit$result$sampler_diagnostics(format = "df")
    summary_diagnostics <- fit$result$diagnostic_summary()
    summary <- fit$result$summary()

    out <- list(
      draws = draws,
      diagnostics = diagnostics,
      summary_diagnostics = summary_diagnostics,
      summary = summary
    )

    # Run diagnostic tests, and message if a flag doesn't pass. Still return
    # the same data
  }

  return(out)
}

#' Get MCMC options
#'
#' @description
#' This function returns a list of MCMC settings to pass to the
#' `cmdstanr::sample()` function to fit the model. The default settings are
#' specified for production-level runs, consider adjusting to optimize
#' for speed while iterating.
#'
#'
#' @param iter_warmup integer indicating the number of warm-up iterations,
#' default is `750`
#' @param iter_sampling integer indicating the number of sampling iterations,
#' default is `500`
#' @param n_chains integer indicating the number of MCMC chains to run, default
#' is `4`
#' @param seed set of integers indicating the random seed of the stan sampler,
#' default is `123`
#' @param adapt_delta float between 0 and 1 indicating the average acceptance
#' probability, default is `0.95`
#' @param max_treedepth integer indicating the maximum tree depth of the
#' sampler, default is 12
#' @param compute_likelihood integer indicating whether or not to compute the
#' likelihood using the data, default is `1` which will fit the model to the
#' data. If set to 0, the model will sample from the prior only
#'
#' @return a list of mcmc settings with the values given by the  function
#' arguments
#' @export
#'
#' @examples
#' mcmc_settings <- get_mcmc_options()
get_mcmc_options <- function(
    iter_warmup = 750,
    iter_sampling = 500,
    n_chains = 4,
    seed = 123,
    adapt_delta = 0.95,
    max_treedepth = 12,
    compute_likelihood = 1) {
  mcmc_settings <- list(
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    n_chains = n_chains,
    seed = seed,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    compute_likelihood = compute_likelihood
  )

  return(mcmc_settings)
}

#' Get model specificaitons
#' @description
#' This function returns a nested list containing the model specifications
#' in the function arguments. All defaults are set for the case of fitting a
#' post-omicron COVID-19 model with joint inference of hospital admissions
#' and data on wastewater viral concentrations
#'
#'
#' @param forecast_date a character string in ISO8 format (YYYY-MM-DD)
#' indicating the date that the forecast is to be made. Default is
#' @param calibration_time integer indicating the number of days to calibrate
#' the model for, default is `90`
#' @param forecast_horizon integer indicating the number of days, including the
#' forecast date, to produce forecasts for, default is `28`
#' @param generation_interval vector of a simplex (must sum to 1) describing
#' the daily probability of onwards transmission, default is package data
#' provided for the COVID-19 generation interval post-Omicron
#' @param inf_to_count_delay vector of a simplex (must sum to 1) describing the
#' daily probability of transitioning from infection to whatever the count
#' variable is, e.g. hospital admissions or cases. Default corresonds to the
#' delay distribution from COVID-19 infection to hospital admission
#' @param infection_feedback_pmf vector of a simplex (must sum to 1) describing
#' the delay from incident infection to feedback in the transmission dynamics.
#' The default is the COVID-19 generation interval
#' @param params a 1 row dataframe of parameters corresponding to model
#' priors and disease/data specific parameters. Default is for COVID-19 hospital
#' admissions and viral concentrations in wastewater
#'
#' @return a list of model specs to be passed to the `get_stan_data()` function
#' @export
#'
#' @examples
#' model_spec_list <- model_spec(forecast_date = "2023-12-06")
model_spec <- function(
    forecast_date,
    calibration_time = 90,
    forecast_horizon = 28,
    generation_interval = wwinference::generation_interval,
    inf_to_count_delay = wwinference::inf_to_hosp,
    infection_feedback_pmf = wwinference::generation_interval,
    params = get_params(
      system.file("extdata", "example_params.toml",
        package = "wwinference"
      )
    )) {
  model_specs <- list(
    forecast_date = forecast_date,
    calibration_time = calibration_time,
    forecast_horizon = forecast_horizon,
    generation_interval = generation_interval,
    inf_to_count_delay = inf_to_hosp,
    infection_feedback_pmf = infection_feedback_pmf,
    params = params
  )
  return(model_specs)
}
