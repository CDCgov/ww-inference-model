#' Example wastewater dataset.
#'
#' A dataset containing the simulated wastewater concentrations
#' (labeled here as `log_genome_copies_per_ml`) by sample collection date
#' (`date`), the site where the sample was collected (`site`) and the lab
#' where the samples were processed (`lab`). Additional columns that are
#' required attributes needed for the model are the limit of detection for
#' that lab on each day (labeled here as `log_lod`) and the population size of
#' the wastewater catchment area represented by the wastewater concentrations
#' in each `site`.
#'
#' This data is generated via the default values in the
#'  `generate_simulated_data()` function. They represent the bare minumum
#'  required fields needed to pass to the model, and we recommend that users
#'  try to format their own data to match this format.
#'
#' The variables are as follows:
#'
#' @format ## ww_data
#' A tibble with 102 rows and 6 columns
#' \describe{
#'   \item{date}{Sample collection date, formatted in ISO8601 standards as
#'   YYYY-MM-DD}
#'   \item{site}{The wastewater treatment plant where the sample was collected}
#'   \item{lab}{The lab where the sample was processed}
#'   \item{log_genome_copies_per_ml}{The natural log of the wastewater
#'   concentration measured on the date specified, collected in the site
#'   specified, and processed in the lab specified. The package expects
#'   this quantity in units of log estimated genome copies per mL.}
#'   \item{log_lod}{The log of the limit of detection in the site and lab on a
#'   particular day of the quantification device (e.g. PCR).  This should be in
#'    units of log estimated genome copies per mL.}
#'   \item{site_pop}{The population size of the wastewater catchment area
#'   represented by the site variable}
#'   \item{location}{ A string indicating the location that all of the
#'   data is coming from. This is not a necessary column, but instead is
#'   included to more realistically mirror a typical workflow}
#'   }
#' @source vignette_data.R
"ww_data"

#' Example evaluation wastewater dataset.
#'
#' A dataset containing the simulated retrospective wastewater concentrations
#' (labeled here as `log_genome_copies_per_ml_eval`) by sample collection date
#' (`date`), the site where the sample was collected (`site`) and the lab
#' where the samples were processed (`lab`). Additional columns that are
#' required attributes needed for the model are the limit of detection for
#' that lab on each day (labeled here as `log_lod`) and the population size of
#' the wastewater catchment area represented by the wastewater concentrations
#' in each `site`.
#'
#' This data is generated via the default values in the
#'  `generate_simulated_data()` function. They represent the bare minumum
#'  required fields needed to pass to the model, and we recommend that users
#'  try to format their own data to match this format.
#'
#' The variables are as follows:
#'
#' @format ## ww_data_eval
#' A tibble with 126 rows and 6 columns
#' \describe{
#'   \item{date}{Sample collection date, formatted in ISO8601 standards as
#'   YYYY-MM-DD}
#'   \item{site}{The wastewater treatment plant where the sample was collected}
#'   \item{lab}{The lab where the sample was processed}
#'   \item{log_genome_copies_per_ml_eval}{The natural log of the wastewater
#'   concentration measured on the date specified, collected in the site
#'   specified, and processed in the lab specified. The package expects
#'   this quantity in units of log estimated genome copies per mL.}
#'   \item{log_lod}{The log of the limit of detection in the site and lab on a
#'   particular day of the quantification device (e.g. PCR).  This should be in
#'    units of log estimated genome copies per mL.}
#'   \item{site_pop}{The population size of the wastewater catchment area
#'   represented by the site variable}
#'   \item{location}{ A string indicating the location that all of the
#'   data is coming from. This is not a necessary column, but instead is
#'   included to more realistically mirror a typical workflow}
#'   }
#' @source vignette_data.R
"ww_data_eval"


#' Example hospital admissions dataset
#'
#'  A dataset containing the simulated daily hospital admissions
#' (labeled here as `daily_hosp_admits`) by date of admission (`date`).
#'  Additional columns that are required are the population size of the
#'  population contributing to the hospital admissions. It is assumed that
#'  the wastewater sites are subsets of this larger population, which
#'  is in the package data assumed to be from a hypothetical US state.
#'  The data generated are daily hospital admissions but they could be any other
#'  epidemiological count dataset e.g. cases. This data should only contain
#'  hospital admissions that would have been available as of the date that
#'  the forecast was made. We recommend that users try to format their data
#'  to match this format.
#'
#' This data is generated via the default values in the
#'  `generate_simulated_data()` function. They represent the bare minimum
#'  required fields needed to pass to the model, and we recommend that users
#'  try to format their own data to match this format.
#'
#' The variables are as follows:
#' \describe{
#'   \item{date}{Date the hospital admissions occurred, formatte din ISO8601
#'   standatds as YYYY-MM-DD}
#'   \item{daily_hosp_admits}{The number of individuals admitted to the
#'   hospital on that date, available as of the forecast date}
#'   \item{state_pop}{The number of people contributing to the daily hospital
#'   admissions}
#'   \item{location}{ A string indicating the location that all of the
#'   data is coming from. This is not a necessary column, but instead is
#'   included to more realistically mirror a typical workflow}
#'   }
#' @source vignette_data.R
"hosp_data"

#' Global R(t) estimate dataset
#'
#'  A dataset containing the the simulated unadjusted daily R(t) estimate
#'  as well as the realized R(t) estimate for the overall "global" population.
#'  This realized R(t) is generated by summing the incident infections from
#'  each of the subpopulations (the populations in each wastewater catchment
#'  area plus an additional one representing the remainder of the population).
#'  The R(t) is then back-calculated by dividing the total incident infection
#'  time seres by the convolution of the infection time series and the
#'  generation interval.
#'
#' This data is generated via the default values in the
#'  `generate_simulated_data()` function.
#'
#' The variables are as follows:
#' \describe{
#'   \item{unadj_rt}{The daily global unadjusted R(t). The log of the
#'   unadjusted R(t) is the central value around which the subpopulation
#'   R(t) values are drawn from each week.}
#'   \item{realized_rt}{The daily value of the realized global R(t). This is
#'   calculated by summing up the subpopulation level incident infections and
#'   then dividing by the convolution of the total incident infections and the
#'   generation interval to get the R(t) value that corresponds to the
#'   global incident infections}
#'   \item{t}{The time index in days}
#'   \item{date}{Date the hospital admissions occurred, formatte din ISO8601
#'   standatds as YYYY-MM-DD}
#'   }
#' @source vignette_data.R
"true_global_rt"


#' Example hospital admissions dataset for evaluation
#'
#'  A dataset containing the simulated daily hospital admissions that the model
#'  will be evaluated against (labeled here as `daily_hosp_admits_for_eval`)
#'  by date of admission (`date`). This data is not needed to fit the model,
#'  but is used in the Getting Started vignette to demonstrate the forecasted
#'  hospital admissions compared to those later observed.
#'
#' This data is generated via the default values in the
#'  `generate_simulated_data()` function.
#'
#' The variables are as follows:
#' \describe{
#'   \item{date}{Date the hospital admissions occurred, formatte din ISO8601
#'   standatds as YYYY-MM-DD}
#'   \item{daily_hosp_admits_for_eval}{The number of individuals admitted to the
#'   hospital on that date, available beyond the forecast date for evaluating
#'   the forecasted hospital admissions}
#'   \item{state_pop}{The number of people contributing to the daily hospital
#'   admissions}
#'   }
#' @source vignette_data.R
"hosp_data_eval"


#' Example subpopulation level hospital admissions dataset
#'
#'  A dataset containing the simulated daily hospital admissions
#' (labeled here as `daily_hosp_admits`) by date of admission (`date`) in
#'  each subpopulation.
#'  Additional columns that are the population size of the
#'  population contributing to the hospital admissions. In this instance,
#'  the subpopulations here are each of the wastewater catchment areas plus
#'  an additional subpopulation for the portion of the population not captured
#'  by wastewater surveillance. The data generated are daily hospital
#'  admissions but they could be any other epidemiological count dataset e.g.
#'  cases. This data should only contain hospital admissions that would have
#'  been available as of the date that the forecast was made.
#'
#' This data is generated via the default values in the
#'  `generate_simulated_data()` function.
#'
#' The variables are as follows:
#' \describe{
#'   \item{date}{Date the hospital admissions occurred, formatted in ISO8601
#'   standards as YYYY-MM-DD}
#'   \item{subpop_name}{A string indicating the subpopulation the hospital
#'   admissiosn corresponds to. This is either a wastewater site, or the
#'   remainder of the population}
#'   \item{daily_hosp_admits}{The number of individuals admitted to the
#'   hospital on that date, available as of the forecast date}
#'   \item{subpop_pop}{The number of people contributing to the daily hospital
#'   admissions in each subpopulation}
#'   }
#' @source vignette_data.R
"subpop_hosp_data"


#' Example subpopulation level retrospective hospital admissions dataset
#'
#'  A dataset containing the simulated daily hospital admissions
#' (labeled here as `daily_hosp_admits`) by date of admission (`date`) in
#'  each subpopulation observed retrospectively.
#'  Additional columns that are required are the population size of the
#'  population contributing to the hospital admissions. In this instance,
#'  the subpopulations here are each of the wastewater catchment areas plus
#'  an additional subpopulation for the portion of the population not captured
#'  by wastewater surveillance. The data generated are daily hospital
#'  admissions but they could be any other epidemiological count dataset e.g.
#'  cases.This data should contain hospital admissions retrospectively beyond
#'  the forecast date in order to evaluate the forecasts.
#'
#'  This data is generated via the default values in the
#'  `generate_simulated_data()` function. They represent the bare minimumum
#'  required fields needed to pass to the model, and we recommend that users
#'  try to format their own data to match this format.
#'
#' The variables are as follows:
#' \describe{
#'   \item{date}{Date the hospital admissions occurred, formatted in ISO8601
#'   standards as YYYY-MM-DD}
#'   \item{subpop_name}{A string indicating the subpopulation the hospital
#'   admissions corresponds to. This is either a wastewater site, or the
#'   remainder of the population}
#'   \item{daily_hosp_admits_for_eval}{The number of individuals admitted to the
#'   hospital on that date, available as of the forecast date}
#'   \item{subpop_pop}{The number of people contributing to the daily hospital
#'   admissions in each subpopulation}
#'   }
#' @source vignette_data.R
"subpop_hosp_data_eval"


#' COVID-19 post-Omicron generation interval probability mass function
#'
#' \describe{
#' A vector that sums to 1, with each element representing the daily
#' probability of secondary onward transmission occurring on that day. The
#' first element of this vector represents the day after primary transmission
#' occurred, it is assumed to be impossible for primary and secondary
#' transmission to occur on the same day.
#' }
#' @source covid_pmfs.R
"default_covid_gi"

#' COVID-19 time delay distribution from infection to hospital admission
#'
#' \describe{
#' A vector that sums to 1, with each element representing the daily
#' probabilty of transitioning from infected to hospitalized, conditioned on
#' being infected and eventually ending up hospitalized. The first element
#' represents the probability of being infected and admitted to the hospital
#' on the same day
#' }
#' @source covid_pmfs.R
"default_covid_inf_to_hosp"
