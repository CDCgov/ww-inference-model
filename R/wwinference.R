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
#' @param forecast_date a character string in ISO8601 format (YYYY-MM-DD)
#' indicating the date that the forecast is to be made. Default is NULL
#' @param calibration_time integer indicating the number of days to calibrate
#' the model for, default is `90`
#' @param forecast_horizon integer indicating the number of days, including the
#' forecast date, to produce forecasts for, default is `28`
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
#' @param dist_matrix Distance matrix, n_sites x n_sites, passed to a
#' distance-based correlation function for epsilon. If NULL, use an independence
#' correlation function, for current implementation (i.e. all sites' epsilon
#' values are independent and identically distributed) .
#' @param corr_structure_switch Integer variable to define the type of
#' correlation matrix structure used.  Input 0 for an iid correlation structure,
#' 1 for an exponential correlation structure based off distance matrix, and 3
#' to use an unstructured, lkj, correlation matrix.
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
                        forecast_date = NULL,
                        calibration_time = 90,
                        forecast_horizon = 28,
                        model_spec = get_model_spec(),
                        mcmc_options = get_mcmc_options(),
                        generate_initial_values = TRUE,
                        compiled_model = compile_model(),
                        dist_matrix = NULL,
                        corr_structure_switch = 0) {
  if (is.null(forecast_date)) {
    cli::cli_abort(
      "The user must specify a forecast date"
    )
  }

  # Check that data is compatible with specifications
  assert_no_dates_after_max(ww_data$date, forecast_date)
  assert_no_dates_after_max(count_data$date, forecast_date)

  input_count_data <- get_input_count_data_for_stan(
    count_data,
    calibration_time
  )
  last_count_data_date <- max(input_count_data$date, na.rm = TRUE)
  first_count_data_date <- min(input_count_data$date, na.rm = TRUE)
  input_ww_data <- get_input_ww_data_for_stan(
    ww_data,
    first_count_data_date,
    last_count_data_date,
    calibration_time
  )
  input_data <- list(
    input_count_data = input_count_data,
    input_ww_data = input_ww_data
  )

  if (corr_structure_switch == 1 && is.null(dist_matrix)) {
    stop(
      "Spatial Components Desired, but Distance Matrix Not Supplied!!!\n
          *distance matrix required for current implementation*"
    )
  }


  # If checks pass, create stan data object
  stan_args <- get_stan_data(
    input_count_data = input_count_data,
    input_ww_data = input_ww_data,
    forecast_date = forecast_date,
    calibration_time = calibration_time,
    forecast_horizon = forecast_horizon,
    generation_interval = model_spec$generation_interval,
    inf_to_count_delay = model_spec$inf_to_count_delay,
    infection_feedback_pmf = model_spec$infection_feedback_pmf,
    params = model_spec$params,
    include_ww = as.numeric(model_spec$include_ww),
    compute_likelihood = as.integer(model_spec$compute_likelihood),
    dist_matrix = dist_matrix,
    corr_structure_switch = corr_structure_switch
  )

  init_lists <- NULL
  if (generate_initial_values) {
    init_lists <- c()
    for (i in 1:mcmc_options$n_chains) {
      init_lists[[i]] <- get_inits_for_one_chain(stan_args, params)
    }
  }


  fit_model <- function(compiled_model,
                        stan_args,
                        model_spec,
                        init_lists) {
    fit <- compiled_model$sample(
      data = stan_args,
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
    stan_args,
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
        to = min(count_data$date) + stan_args$ot + stan_args$ht,
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
        site_pop = stan_args$subpop_size[
          length(unique(stan_args$subpop_size))
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
    max_treedepth = 12) {
  mcmc_settings <- list(
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    n_chains = n_chains,
    seed = seed,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
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
#' @param include_ww Boolean indicating whether or not to include the wastewater
#' data into the model, default is `TRUE` which will get passed to stan and
#' tell the model to evaluate the likelihood with the wastewater data
#' @param compute_likelihood Boolean indicating whether or not to compute the
#' likelihood using the data, default is `TRUE` which will fit the model to the
#' data. If set to `FALSE`, the model will sample from the priors.
#' @param params a list of parameters corresponding to model
#' priors and disease/data specific parameters. Default is for COVID-19 hospital
#' admissions and viral concentrations in wastewater
#'
#' @return a list of model specs to be passed to the `get_stan_data()` function
#' @export
#'
#' @examples
#' model_spec_list <- get_model_spec()
get_model_spec <- function(
    generation_interval = wwinference::default_covid_gi,
    inf_to_count_delay = wwinference::default_covid_inf_to_hosp,
    infection_feedback_pmf = wwinference::default_covid_gi,
    include_ww = TRUE,
    compute_likelihood = TRUE,
    params = get_params(
      system.file("extdata", "example_params.toml",
        package = "wwinference"
      )
    )) {
  model_specs <- list(
    generation_interval = generation_interval,
    inf_to_count_delay = inf_to_count_delay,
    infection_feedback_pmf = infection_feedback_pmf,
    include_ww = include_ww,
    compute_likelihood = compute_likelihood,
    params = params
  )
  return(model_specs)
}
