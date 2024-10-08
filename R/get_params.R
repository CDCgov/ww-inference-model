## Get model parameters (priors and set params of set distributions)------------
#' @title Get parameters for model run
#'
#' @param param_file Path to a `.toml` file defining
#' parameter values.
#'
#' @return a flat list mapping parameter names to corresponding
#' values.
#' @export
#'
#' @examples
#' params <- get_params(param_file = fs::path_package("extdata",
#'   "example_params.toml",
#'   package = "wwinference"
#' ))
get_params <- function(param_file) {
  paramlist <- RcppTOML::parseTOML(param_file)
  validate_paramlist(paramlist)

  flat_paramlist <- c(
    paramlist$continuous_distribution_parameters,
    paramlist$timescale,
    paramlist$infection_process,
    paramlist$hospital_admission_observation_process,
    paramlist$wastewater_observation_process
  )

  return(flat_paramlist)
}

#' Validate a parameter list
#'
#' @param paramlist The parameter list to validate
#'
#' @return the parameter list, on success,
#' or raise an error
#' @export
validate_paramlist <- function(paramlist) {
  expected_sections <- c(
    "continuous_distribution_parameters",
    "hospital_admission_observation_process",
    "infection_process",
    "timescale",
    "wastewater_observation_process"
  )


  missing_sections <- setdiff(
    names(paramlist),
    expected_sections
  )

  if (length(missing_sections) > 0) {
    cli::cli_abort(
      paste0(
        "Parameter list missing expected ",
        "section(s) {missing_sections}"
      )
    )
  }

  ## additional validation logic can go here
  return(paramlist)
}
