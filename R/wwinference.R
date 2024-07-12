#' @title Joint inference of count data (e.g. cases/admissions) and wastewater
#' data
#'
#' @description
#' Provides a user friendly interface around package functionality
#' to produce estimates, nowcasts, and forecasts pertaining to user-specified
#' delay distributions, set parameters, and priors that can be modified to
#' handledifferent types of "global" count data and "local" wastewater
#' concentrationdata using a Bayesian hierarchical framework applied to the two
#' distinctdata sources. By default the model assumes a fixed generation
#' interval and delay from infection to the event that is counted. See the
#' getting started vignette for an example model specifications fitting
#' COVID-19 hospital admissions from a hypothetical state and wasteawter
#' concentration data from multiple sites within that state.
#'
#' @param ww_data A dataframe containing the pre-processed, site-level
#' wastewater concentration data for a model run. The dataframe must contain
#' the following columns: `date`, `site`, `lab`, `genome_copies_per_ml`,
#' `lab_site_index`, `lod`, `below_lod`, `site_pop` `exclude`
#' @param count_data A dataframe containing the pre-procssed, "global" (e.g.
#' state) daily count data, pertaining to the number of events that are being
#' counted on that day, e.g. number of daily cases or daily hospital admissions.
#' Must contain the following columns: `date`, `count` , `total_pop`
#' @param model_spec The model specification parameters as defined using
#' `get_model_spec()`. The default here pertains to the `forecast_date` in the
#' example data provided by the package, but this should be specified by the
#' user based on the date they are producing a forecast
#' @param mcmc_options The MCMC parameters as defined using
#' `get_mcmc_options()`.
#' @param generate_initial_values Boolean indicating whether or not to specify
#' the initialization of the sampler, default is `TRUE`, meaning that
#' initialization lists will be generated and passed as the `init` argument
#' to the model object [`$sample()`][cmdstanr::model-method-sample] call.
#' function
#' @param compiled_model The pre-compiled model as defined using
#' `compile_model()`
#'
#' @return A nested list of the following items, intended to allow the user to
#' quickly and easily plot results from their inference, while also being able
#' to have the full user functionality of running the model themselves in stan
#' by providing the raw model object and diagnostics. If the model runs, this
#' function will return:
#' `draws_df`: A tibble containing the full set of posterior draws of the
#' estimated, nowcasted, and forecasted: counts, site-level wastewater
#' concentrations, "global"(e.g. state) R(t) estimate, and the  "local" (site +
#' the one auxiliary subpopulation) R(t) estimates. In the instance where there
#' are observations, the data will be joined to each draw of the predicted
#' observation to facilitate plotting.
#' `raw_fit_obj`: The CmdStan object that is returned from the call to
#' `cmdstanr::sample()`. Can be used to access draws, summary, diagnostics, etc.
#' `date_time_spine`: Mapping from time in stan to dates
#' `lab_site_spine`: Mapping from lab_site_index in stan to lab and site
#' `subpop_spine`: Mapping from site index in stan to site
#'
#' If the model fails to run, a list containing the follow will be returned:
#' `error`: the error message provided from stan, indicating why the model
#' failed to run. Note, the model might still run and produce draws even if it
#' has major model issues. We recommend the user always run the
#' `check_diagnostics()` function on the `diagnostic_summary` as part of any
#' pipeline to ensure model convergence.
#' @export
#'
wwinference <- function(ww_data,
                        count_data,
                        model_spec = wwinference::get_model_spec(
                          forecast_date =
                            "2023-12-06"
                        ),
                        mcmc_options = wwinference::get_mcmc_options(),
                        generate_initial_values = TRUE,
                        compiled_model = wwinference::compile_model()) {
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

  init_lists <- NULL
  if (generate_initial_values) {
    init_lists <- c()
    for (i in 1:mcmc_options$n_chains) {
      init_lists[[i]] <- get_inits(stan_data, params)
    }
  }


  fit_model <- function(compiled_model,
                        standata,
                        model_spec,
                        init_lists) {
    fit <- compiled_model$sample(
      data = stan_data,
      init = init_lists,
      seed = mcmc_options$seed,
      iter_sampling = mcmc_options$iter_sampling,
      iter_warmup = mcmc_options$iter_warmup,
      max_treedepth = mcmc_options$max_treedepth,
      chains = mcmc_options$n_chains,
      parallel_chains = mcmc_options$n_chains
    )

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
    message(fit$error[[1]])
  } else {
    # This is a bit messy, but get the spines needed to map stan data to
    # the real data
    # Time index to date
    date_time_spine <- tibble::tibble(
      date = seq(
        from = min(count_data$date),
        to = min(count_data$date) + stan_data$ot + stan_data$ht,
        by = "days"
      )
    ) |>
      dplyr::mutate(t = row_number())
    # Lab-site index to corresponding lab, site, and site population size
    lab_site_spine <- ww_data |>
      dplyr::distinct(site, lab, lab_site_index, site_pop)
    # Site index to corresponding site and subpopulation size
    subpop_spine <- ww_data |>
      dplyr::distinct(site, site_index, site_pop) |>
      dplyr::mutate(site = as.factor(site)) |>
      dplyr::bind_rows(tibble::tibble(
        site = "remainder of pop",
        site_index = max(ww_data$site_index) + 1,
        site_pop = stan_data$subpop_size[
          length(unique(stan_data$subpop_size))
        ]
      ))

    draws <- get_draws_df(
      ww_data = ww_data,
      count_data = count_data,
      fit_obj = fit,
      date_time_spine = date_time_spine,
      lab_site_spine = lab_site_spine,
      subpop_spine = subpop_spine
    )
    summary_diagnostics <- fit$result$diagnostic_summary()
    convergence_flag_df <- get_model_diagnostic_flags(
      stan_fit_object =
        fit$result
    )

    out <- list(
      draws_df = draws,
      raw_fit_obj = fit$result,
      date_time_spine = date_time_spine,
      lab_site_spine = lab_site_spine,
      subpop_spine = subpop_spine
    )

    # Message if a flag doesn't pass. Still return
    # the same data, but we want user to know the issue
    if (any(convergence_flag_df[1, ])) {
      warning("Model flagged for convergence issues, run model diagnostics
      on the output stanfit object for more information")
    }
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


#' Get model specifications
#' @description
#' This function returns a nested list containing the model specifications
#' in the function arguments. All defaults are set for the case of fitting a
#' post-omicron COVID-19 model with joint inference of hospital admissions
#' and data on wastewater viral concentrations
#'
#'
#' @param forecast_date a character string in ISO8601 format (YYYY-MM-DD)
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
#' model_spec_list <- get_model_spec(forecast_date = "2023-12-06")
get_model_spec <- function(
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
