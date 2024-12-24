library(scoringutils)
# Using scoringutils 2.0.0

# The sample data already has the prediction as `predictions`, the evaluation
# or "truth" data as `observed` and the samples/draws as `sample_id`. The other
# columns are additional metadata that we can specify later in
# `set_forecast_unit`
sample_data <- scoringutils::example_sample_continuous

# Step 1: convert to quantiles usign scoringutils
quantiled_data <- scoringutils::as_forecast_quantile(
  sample_data,
  # These are the HUb quantiles
  probs = c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99),
  type = 7
)

# Step 2: Create a forecast object from the quantile forecast data.
# Note, here we will replace the `set_forecast_unit` entries with
# the metadata relevant to the forecasts e.g. date, model, forecast_date,
# or generative model, however these are labeled.
forecast_object <- quantiled_data |>
  scoringutils::set_forecast_unit(
    c("location", "target_end_date", "target_type", "horizon", "model")
  )

# Step 3: Score the object  (you actually don't have to do this to get
# the quantile coverage but this will give interval and quantile coverage for
# the 50th and 90th quantiles)
scores <- scoringutils::score(forecast_object)

# Step 4: Compute coverage summarized across all observations
coverage <- scoringutils::get_coverage(forecast_object)

# Step 5 Plot the coverage
coverage |> scoringutils::plot_quantile_coverage()
