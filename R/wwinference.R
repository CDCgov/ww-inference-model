wwinference <- function(ww_data,
                        count_data,
                        model_spec = list(
                          forecast_date = "2023-12-06",
                          calibration_time = 90,
                          forecast_horizon = 28,
                          generation_interval =
                            wwinference::generation_interval,
                          inf_to_count_delay = wwinference::inf_to_hosp,
                          infection_feedback_pmf =
                            wwinference::generation_interval,
                          params = get_params(
                            system.file("extdata", "example_params.toml",
                              package = "wwinference"
                            )
                          ),
                          # Default MCMC settings
                          iter_warmup = 750,
                          iter_sampling = 500,
                          n_chains = 4,
                          seed = 123,
                          adapt_delta = 0.95,
                          max_treedepth = 12,
                          # Default fitting to data
                          compute_likelihood = 1
                        ),
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
