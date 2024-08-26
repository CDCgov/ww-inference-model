# `wwinference`: joint inference and forecasting from wastewater and epidemiological indicators

> [!CAUTION]
> This project is a work-in-progress. Despite this project's early stage, all development is in public as part of the Center for Forecasting and Outbreak Analytics' goals around open development. Questions and suggestions are welcome through GitHub issues or a PR.
>

## Overview

This project is an in-development R package, `{wwinference}` that estimates latent incident infections from wastewater concentration data and data on epidemiological indicators, with an initial assumed structure that the wastewater concentration data comes from subsets of the population contributing to the "global" epidemiological indicator data, such as hospital admissions.
In brief, our model builds upon [EpiNow2](https://github.com/epiforecasts/EpiNow2/tree/main), a widely used [R](https://www.r-project.org/) and [Stan](https://mc-stan.org/) package for Bayesian epidemiological inference.
We modify EpiNow2 to add model for the observed viral RNA concentration in wastewater, adding hierarchical structure to link the subpopulations represented by the osberved wastewater concentrations in each wastewater catchment area.
See our Model Definition page for a mathematical description of the generative model, and the Getting Stated vignette to see an example of how to run the inference model on simulated data.

The intention is for {wwinference} to provide a user-friendly R-package interface for running forecasting models that use wastewater concentrations combined with other more traditional epidemiological signals such as cases or hospital admissions. It aims to be a re-implementation of the modeling components contained in the [wastewater-informed-covid-forecasting](https://github.com/CDCgov/wastewater-informed-covid-forecasting) project repository, with
an emphasis here on making it easier for users to supply their own data.

We recommend reading the [model definition](model_definition.md) to learn more about how the model is structured and running the ["Getting Started" vignette](vignettes/wwinference.Rmd) for an example of how to fit the model to simulated data of COVID-19 hospital admissions and wastewater concentrations.
This will help make clear the data requirements and how to structure this data to fit the model.

## Project Admins
- Kaitlyn Johnson (kaitejohnson)
- Dylan Morris (dylanhmorris)
- Sam Abbott (seabbs)
- Damon Bayer (damonbayer)

# Installing and running code

## Install R
To run our code, you will need a working installation of [R](https://www.r-project.org/) (version `4.3.0` or later). You can find instructions for installing R on the official [R project website](https://www.r-project.org/).

## Install `cmdstanr` and `CmdStan`
We do inference from our models using [`CmdStan`](https://mc-stan.org/users/interfaces/cmdstan) (version `2.35.0` or later) via its R interface [`cmdstanr`](https://mc-stan.org/cmdstanr/) (version `0.8.0` or later).

Open an R session and run the following command to install `cmdstanr` per that package's [official installation guide](https://mc-stan.org/cmdstanr/#installation).

```R
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```

`cmdstanr` provides tools for installing `CmdStan` itself. First check that everything is properly configured by running:

```R
cmdstanr::check_cmdstan_toolchain()
```

You should see the following:
```
The C++ toolchain required for CmdStan is setup properly!
```

If you do, you can then install `CmdStan` by running:
```R
cmdstanr::install_cmdstan()
```
If installation succeeds, you should see a message like the following:
```
CmdStan path set to: {a path on your file system}
```

If you run into trouble, consult the official [`cmdstanr`](https://mc-stan.org/cmdstanr/index.html) website for further installation guides and help.

## Install `wwinference`

Once `cmdstanr` and `CmdStan` are installed, the next step is to download this repository and install the package, `wwinference`. The package provides tools for specifying and running the model, and installs other needed dependencies.

### User install

The package can be installed directly from github by running the following within an R session:
```R
install.packages("remotes")
remotes::install_github("CDCgov/ww-inference-model")
```
## R dependencies
Installing the project package should take care of almost all dependencies installations. Confirm that package installation has succeeded by running the following within an R session:

```R
library(wwinference)
```

## Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](DISCLAIMER.md)
and [Code of Conduct](code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/main/CONTRIBUTING.md), [public domain notices and disclaimers](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md), and [code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).
