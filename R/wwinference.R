#' @title Joint inference of count data (e.g. cases/admissions) and wastewater
#' data
#'
#' @description
#' Provides a user friendly interface around package functionality
#' to produce estimates, nowcasts, and forecasts pertaining to user-specified
#' delay distributions, set parameters, and priors that can be modified to
#' handle different types of "global" count data and "local" wastewater
#' concentration data using a Bayesian hierarchical framework applied to the two
#' distinct data sources. By default the model assumes a fixed generation
#' interval and delay from infection to the event that is counted. See the
#' getting started vignette for an example model specifications fitting
#' COVID-19 hospital admissions from a hypothetical state and wasteawter
#' concentration data from multiple sites within that state.
#'
#' @param ww_data A dataframe containing the pre-processed, site-level
#' wastewater concentration data for a model run. The dataframe must contain
#' the following columns: `date`, `site`, `lab`, `log_genome_copies_per_ml`,
#' `lab_site_index`, `log_lod`, `below_lod`, `site_pop` `exclude`
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
#' @param fit_opts The fit options, which in this case default to the
#' MCMC parameters as defined using `get_mcmc_options()`. This includes
#' the following arguments, which are passed to
#' [`$sample()`][cmdstanr::model-method-sample]:
#' the number of chains, the number of warmup
#' and sampling iterations, the maximum tree depth, the average acceptance
#' probability, and the stan PRNG seed
#' @param generate_initial_values Boolean indicating whether or not to specify
#' the initialization of the sampler, default is `TRUE`, meaning that
#' initialization lists will be generated and passed as the `init` argument
#' to the model object [`$sample()`][cmdstanr::model-method-sample] call.
#' function
#' @param initial_values_seed set of integers indicating the random seed of
#' the R sampler of the initial values, default is `NULL`
#' @param compiled_model The pre-compiled model as defined using
#' `compile_model()`
#'
#' @return An object of the `ww_inference_fit` class containing the following
#' items that are intended to be passed to downstream functions to do things
#' like extract posterior draws, get diangostic behavior, and plot results
#' (for example). If the model runs, this
#' function will return:
#' `fit`: The CmdStan object that is returned from the call to
#' `cmdstanr::sample()`. Can be used to access draws, summary, diagnostics, etc.
#' `raw_input_data`: a list containing the input `ww_data` and the input
#' `count_data` used in the model.
#' `stan_data_list`: a list containing the inputs passed directly to the
#' stan model
#' `fit_opts`: a list of the MCMC specifications passed to stan
#'
#' If the model fails to run, a list containing the follow will be returned:
#' `error`: the error message provided from stan, indicating why the model
#' failed to run. Note, the model might still run and produce draws even if it
#' has major model issues. We recommend the user always run the
#' `check_diagnostics()` function on the `parameter_diagnostics` as part of any
#' pipeline to ensure model convergence.
#' @name wwinference
#' @family diagnostics
#'
#' @export
#' @examples
#' \dontrun{
#' ww_data <- tibble::tibble(
#'   date = rep(seq(
#'     from = lubridate::ymd("2023-08-01"),
#'     to = lubridate::ymd("2023-11-01"),
#'     by = "weeks"
#'   ), 2),
#'   site = c(rep(1, 14), rep(2, 14)),
#'   lab = c(rep(1, 28)),
#'   conc = log(abs(rnorm(28, mean = 500, sd = 50))),
#'   lod = log(c(rep(20, 14), rep(15, 14))),
#'   site_pop = c(rep(2e5, 14), rep(4e5, 14))
#' )
#'
#' ww_data_preprocessed <- preprocess_ww_data(ww_data,
#'   conc_col_name = "conc",
#'   lod_col_name = "lod"
#' )
#' input_ww_data <- indicate_ww_exclusions(ww_data_preprocessed)
#'
#' hosp_data <- tibble::tibble(
#'   date = seq(
#'     from = lubridate::ymd("2023-07-01"),
#'     to = lubridate::ymd("2023-10-30"),
#'     by = "days"
#'   ),
#'   daily_admits = sample(5:70, 122, replace = TRUE),
#'   state_pop = rep(1e6, 122)
#' )
#'
#' input_count_data <- preprocess_count_data(
#'   hosp_data,
#'   "daily_admits",
#'   "state_pop"
#' )
#'
#' generation_interval <- to_simplex(c(0.01, 0.2, 0.3, 0.2, 0.1, 0.1, 0.01))
#' inf_to_count_delay <- to_simplex(c(
#'   rep(0.01, 12), rep(0.2, 4),
#'   rep(0.01, 10)
#' ))
#' infection_feedback_pmf <- generation_interval
#'
#' params <- get_params(
#'   system.file("extdata", "example_params.toml",
#'     package = "wwinference"
#'   )
#' )
#' forecast_date <- "2023-11-06"
#' calibration_time <- 90
#' forecast_horizon <- 28
#' include_ww <- 1
#' ww_fit <- wwinference(input_ww_data,
#'   input_count_data,
#'   model_spec = get_model_spec(
#'     forecast_date = forecast_date,
#'     calibration_time = calibration_time,
#'     forecast_horizon = forecast_horizon,
#'     generation_interval = generation_interval,
#'     inf_to_count_delay = inf_to_coutn_delay,
#'     infection_feedback_pmf = infection_feedback_pmf,
#'     params = params
#'   ),
#'   fit_opts = get_mcmc_options(
#'     iter_warmup = 250,
#'     iter_sampling = 250,
#'     n_chains = 2
#'   )
#' )
#' }
#' @rdname wwinference
#' @aliases wwinference_fit
wwinference <- function(ww_data,
                        count_data,
                        forecast_date = NULL,
                        calibration_time = 90,
                        forecast_horizon = 28,
                        model_spec = get_model_spec(),
                        fit_opts = get_mcmc_options(),
                        generate_initial_values = TRUE,
                        initial_values_seed = NULL,
                        compiled_model = compile_model()) {
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
  raw_input_data <- list(
    input_count_data = input_count_data,
    input_ww_data = input_ww_data
  )

  # If checks pass, create stan data object
  stan_data_list <- get_stan_data(
    input_count_data = input_count_data,
    input_ww_data = input_ww_data,
    forecast_date = forecast_date,
    calibration_time = calibration_time,
    forecast_horizon = forecast_horizon,
    generation_interval = model_spec$generation_interval,
    inf_to_count_delay = model_spec$inf_to_count_delay,
    infection_feedback_pmf = model_spec$infection_feedback_pmf,
    params = model_spec$params,
    include_ww = as.integer(model_spec$include_ww),
    compute_likelihood = as.integer(model_spec$compute_likelihood)
  )

  init_lists <- NULL
  if (is.null(initial_values_seed)) {
    initial_values_seed <- sample.int(.Machine$integer.max, 1L)
  }

  if (generate_initial_values) {
    withr::with_seed(initial_values_seed, {
      init_lists <- lapply(
        1:fit_opts$n_chains,
        \(x) {
          get_inits_for_one_chain(stan_data_list)
        }
      )
    })
  }


  # This returns the cmdstan object if the model runs, and result = NULL if
  # the model errors
  safe_fit_model <- purrr::safely(fit_model)

  fit <- safe_fit_model(
    compiled_model = compiled_model,
    stan_data_list = stan_data_list,
    fit_opts = fit_opts,
    init_lists = init_lists
  )

  if (!is.null(fit$error)) { # If the model errors, return the error message
    return(fit$error)
  } else {
    convergence_flag_df <- get_model_diagnostic_flags(fit$result)

    out <- list(
      fit = fit,
      raw_input_data = raw_input_data,
      stan_data_list = stan_data_list,
      fit_opts = fit_opts
    )

    # Message if a flag doesn't pass. Still return
    # the same data, but we want user to know the issue
    if (any(convergence_flag_df[1, ])) {
      warning("Model flagged for convergence issues, run model diagnostics
      on the output stanfit object for more information")
    }
  }

  # Constructs the wwinference_fit class object
  do.call(new_wwinference_fit, out)
}

#' Constructor for the `wwinference_fit` class
#' @param fit The CmdStan object that is the output of fitting the model
#' @param raw_input_data A list containing all the data passed to stan
#' for fitting the model
#' @param stan_data_list A list containing the inputs passed directly to the
#' stan model
#' @param fit_opts A list of the the fitting options, in this case the
#' MCMC specifications passed to stan
#' @return An object of the `wwinference_fit` class.
#' @noRd
#'
new_wwinference_fit <- function(
    fit,
    raw_input_data,
    stan_data_list,
    fit_opts) {
  # Checking
  stopifnot(
    inherits(fit$result, what = "CmdStanFit"),
    inherits(raw_input_data, "list"),
    inherits(stan_data_list, "list"),
    inherits(fit_opts, "list")
  )

  structure(
    list(
      fit = fit,
      raw_input_data = raw_input_data,
      stan_data_list = stan_data_list,
      fit_opts = fit_opts
    ),
    class = "wwinference_fit"
  )
}

#' @param x,object Object of class `wwinference_fit`
#' @param ... Additional parameters passed to the corresponding method
#' @export
#' @rdname wwinference
#' @return
#' - The print method prints out information about the model and
#' returns the object invisibly.
print.wwinference_fit <- function(x, ...) {
  cat("wwinference_fit object\n")
  cat("N of WW sites              :", x$stan_data_list$n_ww_sites, "\n")
  cat("N of unique lab-site pairs :", x$stan_data_list$n_ww_lab_sites, "\n")
  cat("State population           :", formatC(
    x$stan_data_list$state_pop,
    format = "d"
  ), "\n")
  cat("N of weeks                 :", x$stan_data_list$n_weeks, "\n")
  cat("--------------------\n")
  cat("For more details, you can access the following:\n")
  cat(" - `$fit` for the CmdStan object\n")
  cat(" - `$raw_input_data` for the input data\n")
  cat(" - `$stan_data_list` for the stan data arguments\n")
  cat(" - `$fit_opts` for the fitting options\n")
  invisible(x)
}

#' @export
#' @rdname wwinference
#' @return
#' - The summary method returns the outcome from the
#' `$summary` ([cmdstanr::summary()]) function.
summary.wwinference_fit <- function(object, ...) {
  object$fit$result$summary(...)
}


#' Model fitting function
#' @param compiled_model The compiled model object
#' @param stan_data_list The list of data to pass to stan
#' @param fit_opts The fitting specifications
#' @param init_lists A list of initial values for the sampler
#' @return The fit object from the model
#' @noRd
fit_model <- function(compiled_model,
                      stan_data_list,
                      fit_opts,
                      init_lists) {
  fit <- compiled_model$sample(
    data = stan_data_list,
    init = init_lists,
    seed = fit_opts$seed,
    iter_sampling = fit_opts$iter_sampling,
    iter_warmup = fit_opts$iter_warmup,
    max_treedepth = fit_opts$max_treedepth,
    chains = fit_opts$n_chains,
    parallel_chains = fit_opts$parallel_chains,
    show_messages = fit_opts$show_messages,
    refresh = fit_opts$refresh,
    save_latent_dynamics = fit_opts$save_latent_dynamics,
    output_dir = fit_opts$output_dir,
    output_basename = fit_opts$output_basename,
    sig_figs = fit_opts$sig_figs,
    chain_ids = fit_opts$chain_ids,
    threads_per_chain = fit_opts$threads_per_chain,
    opencl_ids = fit_opts$opencl_ids,
    save_warmup = fit_opts$save_warmup,
    thin = fit_opts$thin,
    adapt_engaged = fit_opts$adapt_engaged,
    step_size = fit_opts$step_size,
    metric = fit_opts$metric,
    metric_file = fit_opts$metric_file,
    inv_metric = fit_opts$inv_metric,
    init_buffer = fit_opts$init_buffer,
    term_buffer = fit_opts$term_buffer,
    window = fit_opts$window,
    fixed_param = fit_opts$fixed_param,
    show_exceptions = fit_opts$show_exceptions,
    diagnostics = fit_opts$diagnostics,
    save_metric = fit_opts$save_metric,
    save_cmdstan_config = fit_opts$save_cmdstan_config
  )

  return(fit)
}


#' Get MCMC options
#'
#' @description
#' This function returns a list of MCMC settings to pass to the
#' `cmdstanr::sample()` function to fit the model. The default settings are
#' specified for production-level runs, consider adjusting to optimize
#' for speed while iterating. All input arguments to `cmdstanr::sample()`
#' are configurable by the user. See `cmdstanr::sample()` documentation
#' https://mc-stan.org/cmdstanr/reference/model-method-sample.html for a full
#' description of each argument.
#'
#'
#' @param iter_warmup integer indicating the number of warm-up iterations,
#' default is `750`
#' @param iter_sampling integer indicating the number of sampling iterations,
#' default is `500`
#' @param n_chains integer indicating the number of MCMC chains to run, default
#' is `4`
#' @param parallel_chains integer indicating the number of chains to run
#' in parallel, default is `4`
#' @param seed set of integers indicating the random seed of the Stan sampler,
#' default is NULL
#' @param adapt_delta float between 0 and 1 indicating the average acceptance
#' probability, default is `0.95`
#' @param max_treedepth integer indicating the maximum tree depth of the
#' sampler, default is 12
#' @param show_messages logical indicating whether to print all output
#' during the execution process, default is `TRUE`
#' @param refresh non-negative integer cmdstanr default, default is `NULL`
#' @param save_latent_dynamics logical cmdstanr default, default is `FALSE`
#' @param output_dir string cmdstanr default, default is
#' `getOption("cmdstanr_output_dir")`
#' @param output_basename string cmdstanr default, default is `NULL`
#' @param sig_figs positive integer cmdstanr default, default is  `NULL`
#' @param chain_ids integer vector cdmstanr default, default is
#' `seq_len(n_chains)`,
#' @param threads_per_chain positive integer cmdstanr default, default
#' is  `NULL`
#' @param opencl_ids integer vector of length 2 cmdstanr default, default
#' is `NULL`
#' @param save_warmup logical cmdstanr default, default is `FALSE`
#' @param thin positive integer cmdstanr default, default is `NULL`
#' @param adapt_engaged logical cmdstanr default, default is `TRUE`
#' @param step_size positive real cmdstanr default, default is `NULL`
#' @param metric string cmdstanr default, default is `NULL`
#' @param metric_file character vector cmdstanr default, default is `NULL`
#' @param inv_metric vector, matrix cmdstanr default, default is `NULL``
#' @param init_buffer nonnegative integer cmdstanr default, default is `NULL`
#' @param term_buffer nonnegative integer cmdstanr default, default is `NULL`
#' @param window nonnegative integer cmdstanr default, default is `NULL`
#' @param fixed_param logical cmdstanr default, default is `FALSE`
#' @param show_exceptions logical cmdstanr default, default is  TRUE,
#' @param diagnostics character vector cmdstanr default, default is
#' `c("divergences", "treedepth", "ebfmi")`
#' @param save_metric logical cmdstanr default, default is `NULL`
#' @param save_cmdstan_config logical cmdstanr default, default is `NULL`
#'
#' @return a list of MCMC settings with the values given by the  function
#' arguments
#' @export
#'
#' @examples
#' mcmc_settings <- get_mcmc_options()
get_mcmc_options <- function(
    iter_warmup = 750,
    iter_sampling = 500,
    n_chains = 4,
    parallel_chains = 4,
    seed = NULL,
    adapt_delta = 0.95,
    max_treedepth = 12,
    show_messages = TRUE,
    # CmdstanR default configurations
    refresh = NULL,
    save_latent_dynamics = FALSE,
    output_dir = getOption("cmdstanr_output_dir"),
    output_basename = NULL,
    sig_figs = NULL,
    chain_ids = seq_len(n_chains),
    threads_per_chain = NULL,
    opencl_ids = NULL,
    save_warmup = FALSE,
    thin = NULL,
    adapt_engaged = TRUE,
    step_size = NULL,
    metric = NULL,
    metric_file = NULL,
    inv_metric = NULL,
    init_buffer = NULL,
    term_buffer = NULL,
    window = NULL,
    fixed_param = FALSE,
    show_exceptions = TRUE,
    diagnostics = c("divergences", "treedepth", "ebfmi"),
    save_metric = NULL,
    save_cmdstan_config = NULL) {
  mcmc_settings <- list(
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    n_chains = n_chains,
    seed = seed,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    show_messages = show_messages,
    refresh = refresh,
    save_latent_dynamics = save_latent_dynamics,
    output_dir = output_dir,
    output_basename = output_basename,
    sig_figs = sig_figs,
    parallel_chains = parallel_chains,
    chain_ids = chain_ids,
    threads_per_chain = threads_per_chain,
    opencl_ids = opencl_ids,
    save_warmup = save_warmup,
    thin = thin,
    adapt_engaged = adapt_engaged,
    step_size = step_size,
    metric = metric,
    metric_file = metric_file,
    inv_metric = inv_metric,
    init_buffer = init_buffer,
    term_buffer = term_buffer,
    window = window,
    fixed_param = fixed_param,
    show_exceptions = show_exceptions,
    diagnostics = diagnostics,
    save_metric = save_metric,
    save_cmdstan_config = save_cmdstan_config
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
