#' @keywords internal
"_PACKAGE"

#' @importFrom lubridate ymd
#' @importFrom tidybayes spread_draws stat_halfeye stat_slab
#' @importFrom dplyr filter left_join select pull distinct mutate as_tibble
#' rename ungroup arrange row_number group_by lead lag
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom ggplot2 ggplot facet_wrap geom_line geom_hline geom_point
#' geom_bar theme scale_y_continuous scale_colour_discrete scale_fill_discrete
#' geom_ribbon scale_x_date facet_grid geom_vline labs aes ggtitle
#' @importFrom cmdstanr cmdstan_model
#' @importFrom posterior subset_draws as_draws_list
#' @importFrom fs path_package
#' @importFrom rlang sym
#' @importFrom stats dnbinom dweibull ecdf plogis qlogis rlnorm rnbinom rnorm
#' rt sd time
NULL
