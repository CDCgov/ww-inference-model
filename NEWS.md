# wwinference 0.1.1.99 (dev)
- Fixes the way the wastewater calibration and evaluation datasets are created
to ensure that the dates align properly. ([#256](https://github.com/CDCgov/ww-inference-model/issues/256))

# wwinference 0.1.1
This release includes a change to the default model priors. The previous release had an informative prior on a high magnitude of infection feedback. This release reduces the prior mode and also decreases the degree of prior certainty. This release also includes minor changes to plotting and pre-processing functions designed to make outputs more comprehensive and interpretable.

## User-visible changes
- Add wastewater data into the forecast period to output in `generate_simulated_data()` function and as package data. Also adds subpopulation-level
hospital admissions to output of function and package data. ([#184](https://github.com/CDCgov/ww-inference-model/issues/184))
- `wwinference` now checks whether `site_pop` is fixed per site (see issue [#223](https://github.com/CDCgov/ww-inference-model/issues/226) reported by [@akeyel](https://github.com/akeyel)).
- Add package workflow diagram to readme ([#248](https://github.com/CDCgov/ww-inference-model/issues/248))
- `get_plot_subpop_rt()` now uses a shared y-axis to facilitate comparison of R(t) estimates) ([#245](https://github.com/CDCgov/ww-inference-model/issues/245))

## Internal changes
- Modified the priors on the infection feedback term and the step size of the weekly random walk in the effective reproductive number (issue [#227](https://github.com/CDCgov/ww-inference-model/issues/227)), based on benchmarking results from the evaluation pipeline described in the [PR](https://github.com/CDCgov/ww-inference-model/pull/236) corresponding to this change.
- Updated the workflow for posting the pages artifact to PRs (issue [#229](https://github.com/CDCgov/ww-inference-model/issues/229)(https://github.com/CDCgov/ww-inference-model/issues/229)).
- Modify `plot_forecasted_counts()` so that it does not require an evaluation dataset ([#218](https://github.com/CDCgov/ww-inference-model/pull/218))

# wwinference 0.1.0

This is the first major release, focused on providing an initial version of the package.
Note the package is still flagged as in development, though the authors plan on using it for production work in the coming weeks.
As it's written, the package is intended to allow users to do the following:

- Provide basic functionality to fit the wastewater-informed model to an example fitting COVID-19 hospital admissions and wastewater from a few sites ([#5](https://github.com/CDCgov/ww-inference-model/issues/5))
- Performs basic post-processing and plotting of data and modeled outputs, including calibrated, nowcasted, and forecasted count data (in the example, hospital admissions), wastewater concentrations, global R(t) estimates and subpopulation-level R(t) estimates
- Provide an example in the vignette to fit the model to only the hospital admissions ([#24](https://github.com/CDCgov/ww-inference-model/issues/24))
- Validate input data validation with informative error messaging ([#37](https://github.com/CDCgov/ww-inference-model/issues/37), [#54](https://github.com/CDCgov/ww-inference-model/issues/54))
- Provide a wrapper function to generate forward simulated data with user-specified variables. It calls a number of functions to perform specific model components ([#27](https://github.com/CDCgov/ww-inference-model/issues/27))
- Contains S3 class methods applied to the output of the main model wrapper function, the `wwinference_fit` class object ([#58](https://github.com/CDCgov/ww-inference-model/issues/58)).
- Wastewater concentration data is expected to be in log scale ([#122](https://github.com/CDCgov/ww-inference-model/pull/122)).
