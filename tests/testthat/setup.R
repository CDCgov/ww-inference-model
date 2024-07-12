testthat_stan_include <- function() {
  system.file(
    "stan",
    package = "wwinference"
  )
}

model_file_path_id <- system.file(
  "stan", "wwinference.stan",
  package = "wwinference"
)

cat("\nsetup.R is compiling the stan model in preparation for testing.\n")

# precompiled site-level infection dynamics model
compiled_site_inf_model <- cmdstanr::cmdstan_model(
  model_file_path_id,
  force_recompile = TRUE,
  compile_standalone = TRUE,
  include = testthat_stan_include(),
  dir = tempdir()
)

params <- wwinference::get_params(
  system.file("extdata", "example_params.toml",
    package = "wwinference"
  )
)
