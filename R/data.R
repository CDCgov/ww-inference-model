#' Example wastewater dataset.
#'
#' A dataset containing the simulated wastewater concentrations
#' (labeled here as `genome_copies_per_ml`) by sample collection date (`date`),
#' the site where the sample was collected (`site`) and the lab where the
#' samples were processed (`lab`). Additional columns that are required
#' attributes needed for the model are the limit of detection for that lab on
#' each day (labeled here as `lod`) and the population size of the wastewater
#' catchment area represented by the wastewater concentrations in each `site`.
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
#'   \item{genome_copies_per_ml}{The wastewater concentration measured on the
#'   date specified, collected in the site specified, and processed in the lab
#'   specified. The default parameters assume that this quantity is reported
#'   as the genome copies per mL, on a natural scale.}
#'   \item{lod}{The limit of detection in the site and lab on a particular day
#'   of the quantification device (e.g. PCR). This is also by default reported
#'   in terms of the genome copies per mL.}
#'   \item{site_pop}{The population size of the wastewater catchment area
#'   represented by the site variable}
#'   }
#' @source vignette_data.R
"ww_data"




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
#'  `generate_simulated_data()` function. They represent the bare minumum
#'  required fields needed to pass to the model, and we recommend that users
#'  try to format their own data to match this formate.
#'
#' The variables are as follows:
#' \describe{
#'   \item{date}{Date the hospital admissions occurred, formatte din ISO8601
#'   standatds as YYYY-MM-DD}
#'   \item{daily_hosp_admits}{The number of individuals admitted to the
#'   hospital on that date, available as of the forecast date}
#'   \item{state_pop}{The number of people contributing to the daily hospital
#'   admissions}
#'   }
#' @source vignette_data.R
"hosp_data"

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
"generation_interval"

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
"inf_to_hosp"
