#' Add probability mass functions
#'
#' This function allows the addition of probability mass functions (PMFs) to
#' produce a new PMF. This is useful for example in the context of reporting
#' delays where the PMF of the sum of two Poisson distributions is the
#' convolution of the PMFs.
#'
#' This code was adapted from code written
#' (under an MIT license) as part of the `epinowcast`
#' package (https://github.com/epinowcast/epinowcast)
#'
#' @param pmfs A list of vectors describing the probability mass functions to
#'
#' @return A vector describing the probability mass function of the sum of the
#'
#' @export
#' @examples
#' # Sample and analytical PMFs for two Poisson distributions
#' x <- rpois(10000, 5)
#' xpmf <- dpois(0:20, 5)
#' y <- rpois(10000, 7)
#' ypmf <- dpois(0:20, 7)
#' # Add sampled Poisson distributions up to get combined distribution
#' z <- x + y
#' # Analytical convolution of PMFs
#' conv_pmf <- add_pmfs(list(xpmf, ypmf))
#' conv_cdf <- cumsum(conv_pmf)
#' # Empirical convolution of PMFs
#' cdf <- ecdf(z)(0:42)
#' # Compare sampled and analytical CDFs
#' plot(conv_cdf)
#' lines(cdf, col = "black")
add_pmfs <- function(pmfs) {
  d <- length(pmfs)
  if (d == 1) {
    return(pmfs[[1]])
  }
  if (!is.list(pmfs)) {
    return(pmfs)
  }
  # P(Z = z) = sum_over_x(P(X = x) * P(Y = z - x)) # nolint
  return(
    Reduce(x = pmfs, f = function(conv, pmf) {
      lc <- length(conv)
      wd <- seq_len(lc) - 1
      proc <- numeric(lc + length(pmf))
      for (j in seq_along(pmf)) {
        proc[j + wd] <- proc[j + wd] + pmf[j] * conv
      }
      return(proc)
    })
  )
}

#' @title Get index matrix
#' @description Get a matrix to broadcast a vector from weekly to daily
#' @param n_days number of days we will expand to
#' @param n_weeks number of weeks those days correspond to
#'
#' @return a n_day x n_week matrix for multiplying by weekly estimated
#' value to broadcast it to daily
#' @export
#'
#' @examples
#' ind_m <- get_ind_m(14, 2)
get_ind_m <- function(n_days, n_weeks) {
  ind_m <- matrix(nrow = n_days, ncol = n_weeks)
  for (i in 1:n_days) {
    for (j in 1:n_weeks) {
      if (((i - 1) %/% 7) + 1 == j) {
        ind_m[i, j] <- 1
      } else {
        ind_m[i, j] <- 0
      }
    }
  }

  return(ind_m)
}

#' @title Create a new directory if one doesn't exist
#' @description
#' Function to create a directory for the specified output file path if needed.
#' Does nothing if the target directory already exists.
#'
#'
#' @param output_file_path file path that may or may not need to be created
#'
#' @export
create_dir <- function(output_file_path) {
  if (!file.exists(output_file_path)) {
    fs::dir_create(output_file_path, recurse = TRUE, mode = "0777")
    Sys.chmod(output_file_path, mode = "0777", use_umask = FALSE)
  }
}


#' @title Get the mean of a Normal distribution for a random variable Y
#' needed to ensure that the distribution of X = exp(Y) (which is Log-Normal)
#' has a specified mean and sd.
#' @description
#'  see arithmetic moments here
#' https://en.wikipedia.org/wiki/Log-normal_distribution
#'
#' @param mean target mean for the Log-Normal distribution of X
#' @param sd target sd for the Log-Normal distribution X
#'
#' @return corresponding mean for the underlying Normal
#' distribution of Y = log(X).

#' @export
convert_to_logmean <- function(mean, sd) {
  logmean <- log(mean^2 / sqrt(sd^2 + mean^2))
  return(logmean)
}


#' @title Get the sd of a Normal distribution for a random variable Y
#' needed to ensure that the distribution of X = exp(Y) (which is Log-Normal)
#' has a specified mean and sd.
#' @description see arithmetic moments here
#' https://en.wikipedia.org/wiki/Log-normal_distribution
#'
#' @param mean target mean for the Log-Normal distribution of X
#' @param sd target sd for the Log-Normal distribution of X
#'
#' @return corresponding sd for the underlying Normal distribution of Y = log(X)
#' @export
convert_to_logsd <- function(mean, sd) {
  logsd <- sqrt(log(1 + (sd^2 / mean^2)))
  return(logsd)
}

#' @title Normalize vector to a simplex
#'
#' @param vector numeric vector
#'
#' @return vector whose entries sum to 1
#' @export
#' @examples
#' to_simplex(c(1, 1, 1))
to_simplex <- function(vector) {
  if (any(vector < 0)) {
    cli::cli_abort(
      c(
        "Cannot normalize a vector with ",
        "negative entries to a simplex. ",
        "Got {vector}"
      )
    )
  }
  return(vector / sum(vector))
}

#' Escape brackets returned in a string for passing to glue
#'
#' @param string A string vector containing `{}`
#'
#' @return A string vector where all single brackets are replaced with double
#' brackets
autoescape_brackets <- function(string) {
  return(gsub("\\{|\\}", "", string))
}
