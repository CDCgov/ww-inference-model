#' Compile a stan model while pointing at the package default
#' include directory (`stan`) for #include statements
#'
#' @description
#' This function reads in and optionally compiles a Stan model.
#' It is written to search the installed package `stan` directory
#' for additional source files to include. Within each stan file,
#' use `#include {path to your file with the stan directory}.stan`
#'
#'
#' @details The code for this function has been adapted
#' from code written (under an MIT license) as part of
#' the [`epinowcast`](https://github.com/epinowcast/epinowcast)
#' R package.
#'
#' @param model_filepath path to .stan file defining the model, default is
#' `system.file("stan", "wwinference.stan", package = "wwinference")
#' @param include_paths path(s) to directories to search for files
#' specified in #include statements. Passed to [cmdstanr::cmdstan_model()].
#' Defaults to the `stan` subdirectory of the installed
#' `wwinference` package.
#' @param threads Number of threads to use in model compilation,
#' as an integer. Passed to [cmdstanr::cmdstan_model()].
#' Default `FALSE` (use single-threaded compilation).
#' @param target_dir Directory in which to save the compiled
#' stan model binary. Passed as the `dir` keyword argument to
#' [cmdstanr::cmdstan_model()]. Defaults to a temporary directory
#' for the R session (the output of [tempdir()]).
#' @param stanc_options Options for the stan compiler passed to
#' [cmdstanr::cmdstan_model()], as a list. See that function's
#' documentation for more details. Default `list()` (use default
#' options).
#' @param cpp_options Options for the C++ compiler passed to
#' [cmdstanr::cmdstan_model()], as a list. See that function's
#' documentation for more details. Default `list()` (use default
#' options).
#' @param verbose Write detailed output to the terminal while
#' executing the function? Boolean, default `TRUE`.
#' @param ... Additional keyword arguments passed to
#' [cmdstanr::cmdstan_model()].
#'
#' @return The resulting `cmdstanr` model object, as the output
#' of [cmdstanr::cmdstan_model()].
#' @export
compile_model <- function(model_filepath = system.file("stan",
                            "wwinference.stan",
                            package = "wwinference"
                          ),
                          include_paths = system.file(
                            "stan",
                            package = "wwinference"
                          ),
                          threads = FALSE,
                          target_dir = tempdir(),
                          stanc_options = list(),
                          cpp_options = list(),
                          verbose = TRUE,
                          ...) {
  if (verbose) {
    cli::cli_inform(
      glue::glue(paste0(
        "Using model source file: ",
        "{model_filepath}"
      ))
    )
    cli::cli_inform(
      sprintf(
        "Using include paths: %s",
        toString(include_paths)
      )
    )
  }

  create_dir(target_dir)

  model <- cmdstanr::cmdstan_model(
    model_filepath,
    include_paths = include_paths,
    compile = TRUE,
    stanc_options = stanc_options,
    cpp_options = cpp_options,
    threads = threads,
    dir = target_dir,
    ...
  )

  if (verbose) {
    cli::cli_inform(paste0(
      "Model compiled or loaded successfully; ",
      "model executable binary located at: ",
      "{model$exe_file()}"
    ))
  }

  return(model)
}
