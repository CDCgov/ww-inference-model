#' Get plot of fit and forecasted counts
#'
#' @param draws A dataframe containing the posterior draws with the data joined
#' to it. This is the `draws_df` output of a call to [wwinference()]
#' @param count_data_eval A dataframe containing the count data we will
#' evaluate the forecasts against. Must contain the columns `date` and
#' a column indicating the count data to evaluate against, with the name
#' of that column specified as the `count_data_eval_col_name`
#' @param count_data_eval_col_name string indicating the name of the count
#' data to evaluate against the forecasted count data
#' @param forecast_date A string indicating the date we made the forecast, for
#' plotting, in ISO8601 format YYYY-MM-DD
#' @param count_type A string indicating what data the counts refer to,
#' default is `hospital admissions`
#'
#' @return A ggplot object containing the posterior draw of the estimated,
#' nowcasted, and forecasted counts alongside the data used to
#' calibrate the model and subsequently observed counts (if any) against which
#' to evaluate the forecast performance.

#' @export
#'
get_plot_forecasted_counts <- function(draws,
                                       count_data_eval,
                                       count_data_eval_col_name,
                                       forecast_date,
                                       count_type = "hospital admissions") {
  p <- ggplot(draws |> dplyr::filter(
    name == "pred_counts",
    draw %in% sampled_draws
  )) +
    geom_line(aes(x = date, y = pred_value, group = draw),
      color = "red4", alpha = 0.1, size = 0.2
    ) +
    geom_point(
      data = count_data_eval,
      aes(x = date, y = .data[[count_data_eval_col_name]]),
      shape = 21, color = "black", fill = "white"
    ) +
    geom_point(aes(x = date, y = observed_value)) +
    geom_vline(aes(xintercept = lubridate::ymd(forecast_date)),
      linetype = "dashed"
    ) +
    xlab("") +
    ylab(glue::glue("Daily {count_type}")) +
    ggtitle(glue::glue("Fit and forecasted {count_type}")) +
    scale_x_date(
      date_breaks = "2 weeks",
      labels = scales::date_format("%Y-%m-%d")
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        size = 8, vjust = 1,
        hjust = 1, angle = 45
      ),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(
        size = 10,
        vjust = 0.5, hjust = 0.5
      )
    )
  return(p)
}

#' Get plot of fit and forecasted wastewater concentrations
#'
#' @param draws A dataframe containing the posterior draws with the data joined
#' to it. This is the `draws_df` output of a call to [wwinference()]
#' @param forecast_date A string indicating the date we made the forecast, for
#' plotting, in ISO8601 format YYYY-MM-DD
#'
#' @return a ggplot object containing faceted plots of the wastewaster
#' concentrations in each site and lab combination
#' @export
#'
get_plot_ww_conc <- function(draws,
                             forecast_date) {
  p <- ggplot(draws |> dplyr::filter(
    name == "pred_ww",
    draw %in% sampled_draws
  ) |>
    dplyr::mutate(
      site_lab_name = glue::glue("{subpop}, Lab: {lab}")
    )) +
    geom_line(
      aes(
        x = date, y = log(pred_value),
        color = subpop,
        group = draw
      ),
      alpha = 0.1, size = 0.2,
      show.legend = FALSE
    ) +
    geom_point(aes(x = date, y = log(observed_value)),
      color = "black", show.legend = FALSE
    ) +
    facet_wrap(~site_lab_name, scales = "free") +
    geom_vline(aes(xintercept = lubridate::ymd(forecast_date)),
      linetype = "dashed"
    ) +
    xlab("") +
    ylab("Log(Genome copies/mL)") +
    ggtitle("Lab-site level wastewater concentration") +
    scale_x_date(
      date_breaks = "2 weeks",
      labels = scales::date_format("%Y-%m-%d")
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        size = 8, vjust = 1,
        hjust = 1, angle = 45
      ),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(
        size = 10,
        vjust = 0.5, hjust = 0.5
      )
    )
  return(p)
}

#' Get plot of fit, nowcasted, and forecasted "global" R(t)
#'
#' @param draws A dataframe containing the posterior draws with the data joined
#' to it. This is the `draws_df` output of a call to [wwinference()]
#' @param forecast_date A string indicating the date we made the forecast, for
#' plotting, in ISO8601 format YYYY-MM-DD
#'
#' @return A ggplot object containing the posterior draws of the global R(t)
#' estimate
#' @export
#'
get_plot_global_rt <- function(draws,
                               forecast_date) {
  # R(t) of the hypothetical state
  p <- ggplot(draws |> dplyr::filter(
    name == "global R(t)",
    draw %in% sampled_draws
  )) +
    geom_line(aes(x = date, y = pred_value, group = draw),
      color = "blue4", alpha = 0.1, size = 0.2
    ) +
    geom_vline(aes(xintercept = lubridate::ymd(forecast_date)),
      linetype = "dashed"
    ) +
    geom_hline(aes(yintercept = 1), linetype = "dashed") +
    xlab("") +
    ylab("Global R(t)") +
    ggtitle("Global R(t) estimate") +
    scale_x_date(
      date_breaks = "2 weeks",
      labels = scales::date_format("%Y-%m-%d")
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        size = 8, vjust = 1,
        hjust = 1, angle = 45
      ),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(
        size = 10,
        vjust = 0.5, hjust = 0.5
      )
    )
  return(p)
}

#' Get plot of fit, nowcasted, and forecasted R(t) in each subpopulation
#'
#' @param draws A dataframe containing the posterior draws with the data joined
#' to it. This is `draws_df` output of the call to `wwinference()`
#' @param forecast_date A string indicating the date we made the forecast, for
#' plotting, in ISO8601 format YYYY-MM-DD
#'
#' @return A ggplot object containing faceted plots of the R(t) estimate in each
#' subpopulation (so wastewater sites + those not on wastewater)
#' @export
#'
get_plot_subpop_rt <- function(draws,
                               forecast_date) {
  p <- ggplot(draws |> dplyr::filter(
    name == "subpop R(t)",
    draw %in% sampled_draws
  )) +
    geom_line(
      aes(
        x = date, y = pred_value, group = draw,
        color = subpop
      ),
      alpha = 0.1, size = 0.2,
      show.legend = FALSE
    ) +
    geom_vline(aes(xintercept = lubridate::ymd(forecast_date)),
      linetype = "dashed"
    ) +
    facet_wrap(~subpop, scales = "free") +
    geom_hline(aes(yintercept = 1), linetype = "dashed") +
    xlab("") +
    ylab("Subpopulation R(t)") +
    ggtitle("R(t) estimate of each subpopulation") +
    scale_x_date(
      date_breaks = "2 weeks",
      labels = scales::date_format("%Y-%m-%d")
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        size = 8, vjust = 1,
        hjust = 1, angle = 45
      ),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(
        size = 10,
        vjust = 0.5, hjust = 0.5
      )
    )

  return(p)
}
