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
#' @return vector whos sum adds to 1
#' @export
#' @examples
#' to_simplex(c(1, 1, 1))
#' @noRd
to_simplex <- function(vector) {
  return(vector / sum(vector))
}

#' Drop the first element of a simplex and re-normalize the result to sum to 1.
#'
#' When this vector corresponds to the generation interval distribution, we
#' want to drop this first bin. The renewal equation assumes that same-day
#' infection and onward transmission does not occur, and we assume
#' everything is 1 indexed not 0 indeced. We need to
#' manually drop the first element from the PMF vector.
#'
#' @param x A numeric vector, sums to 1. Corresponds to a discretized PDF or PMF
#'   (usually the GI distribution).
#'
#' @return A numeric vector, sums to 1.
#' @export
#' @examples
#' pmf_orig <- c(0.1, 0.1, 0.1, 0.7)
#' pmf_trunc <- drop_first_and_renormalize(pmf_orig)
drop_first_and_renormalize <- function(x) {
  # Check input sums to 1
  stopifnot(abs(sum(x) - 1) < 1e-8)
  # Drop and renormalize
  y <- x[2:length(x)] / sum(x[2:length(x)])
  vec_outside_tol <- abs(sum(y) - 1L) > 1e-10
  # Normalize until within tolerance
  while (vec_outside_tol) {
    y <- y / sum(y)
    vec_outside_tol <- abs(sum(y) - 1L) > 1e-10
  }
  return(y)
}
