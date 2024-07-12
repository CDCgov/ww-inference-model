## Model overview: wastewater-informed, site-level infection dynamics

This document details the generative model used by `wwinference` to infer global and local infection dynamics from count data (e.g. cases or hospital admissions) and wastewater concentration data, and use that to produce nowcasts and forecasts.
The `wwinference` model assumes that wastewater concentration data is available for one or more "local" sites that represent a subset of the total population that produce the "global" epidemiological indicators (cases or admissions).
In the future, we plan to provide the user with functionality for other types of data structures, e.g. multiple data streams of hospital admissions in addition to multiple wastewater concentration data streams, but for now this structure is the only option provided to the user.
The model is structured into four generative components: [infections](#infection-component), [hierarchical subpopulation-level infection dynamics](#hierarchical-subpopulation-infection-dynamics-component), [hospital admissions](#hospital-admissions-component), and [viral genome concentrations in wastewater](#viral-genome-concentration-in-wastewater-component).
The model imposes a hierarchical structure where the infection dynamics in the subpopulations represented by the wastewater concentration data are assumed to be localized outbreaks similar to one another and centered around the "global" infection dynamics that give rise to the hospital admissions. Note, we will describe the model in terms of the generation of hospital admissions, but the user can choose to replace this with any "count" dataset with a delay distribution from infection to the generation of that count data, e.g. cases would also work well here.

### Model components

Our models are constructed from a set of generative components. These are:

- [**Infection component:**](#infection-component) A renewal model for the infection dynamics, which generates estimates of incident latent infections per capita.
- [**Hierarchical subpopulation-level infection dynamics:**](#hierarchical-subpopulation-infection-dynamics-component) A model for the relation between the subpopulation infection dynamics and the "global" infection dynamics
- [**Hospital admissions component:**](#hospital-admissions-component) A model for the expected number of hospital admissions given incident latent infections per capita.
- [**Viral genome concentration in wastewater:**](#viral-genome-concentration-in-wastewater-component) A model for the expected genome concentration given incident infections per capita.

See the [notation](#appendix-notation) section for an overview of the mathematical notation we use to describe the model components, including how probability distributions are parameterized.

### Infection component

#### Renewal process for incident infections

This component assumes that latent (unobserved) _expected_ incident infections per capita $I(t)$ are generated from a renewal [^Cori][^EpiNow2][^Epidemia] process described by:
$$I(t) = \mathcal{R}(t) \sum_{\tau = 1}^{T_g} I(t-\tau) g(\tau)$$
Where $g(\tau)$ is the discrete generation interval, which describes the distribution of times from incident infection to secondary infection (i.e. infectiousness profile) and $\mathcal{R}(t)$ is the instantaneous reproduction number, representing the expected number of secondary infections occurring at time $t$, divided by the number of currently infected individuals, each scaled by their relative infectiousness at time $t$ [^Gostic2020]. $T_g$ is the maximum generation interval, which is the maximum time from infection to secondary infection that we consider, and is set to 15 days.

This process is initialized by estimating an initial exponential growth[^EpiNow2] of infections for 50 days prior to the calibration start time $t_0$:

$$ I(t) = I_0\exp(rt) $$

where $I_0$ is the initial per capita infection incident infections and $r$ is the exponential growth rate.

#### Instantaneous reproduction number

We decompose the instantaneous reproduction number $\mathcal{R}(t)$ into two components: an _unadjusted_ instantaneous reproduction number $\mathcal{R}^\mathrm{u}(t)$ and a damping term that accounts for the effect of recent infections on the instantaneous reproduction number[^Asher2018].

We assume that the unadjusted reproduction number $\mathcal{R}^\mathrm{u}(t)$ is a piecewise-constant function with weekly change points (i.e., if $t$ and $t'$ are days in the same week, then $\mathcal{R}^\mathrm{u}(t)$ = $\mathcal{R}^\mathrm{u}(t')$ ). To account for the dependence of the unadjusted reproduction number in a given week on the previous week, we use a differenced auto-regressive process for the log-scale reproduction number. A log-scale representation is used to ensure that the reproduction number is positive and so that week-to-week changes are multiplicative rather than additive.

$$
\log[\mathcal{R}^\mathrm{u}(t_i)] \sim \mathrm{Normal}\left(\log[\mathcal{R}^\mathrm{u}(t_{i-1})] + \beta \left(\log[\mathcal{R}^\mathrm{u}(t_{i-2})] - \log[\mathcal{R}^\mathrm{u}(t_1)]\right), \sigma_r \right)
$$

where $t_i$, $t_{i-1}$, and $t_{i-2}$ are days in three successive weeks, $\beta$ is an autoregression coefficient which serves to make week-to-week changes correlated, and $\sigma_r$ determines the overall variation in week-to-week changes.
We bound $\beta$ to be between 0 and 1 so that any changes in trend in $\mathcal{R}^\mathrm{u}(t)$ are damped over time to a degree determined by $\beta$.

The damping term we use is based on Asher et al. 2018[^Asher2018] but extended to be applicable to a renewal process. It assumes that the instantaneous reproduction number is damped by recent infections weighted by the generation interval. This is a simple way to account for the fact that the instantaneous reproduction number is likely to decrease when there are many infections in the population, due to factors such as immunity, behavioral changes, and public health interventions. The damping term is defined as:

$$ \mathcal{R}(t) = \mathcal{R}^\mathrm{u}(t) \exp \left( -\gamma \sum_{\tau = 1}^{T_f}I(t-\tau)g(\tau) \right) $$

where $\gamma$ is the _infection feedback term_ controlling the strength of the damping on $\mathcal{R}(t)$, and the summation is analogous to the "force of infection." See [Prior Distributions](#prior-distributions) below for description of prior choice on $\gamma$.

### Hierarchical subpopulation infection dynamics component

The structure of this model assumes that we have hospital admissions data coming from a larger "global" population (e.g. an entire state or county) and localized wastewater concentration measurements coming from subsets of the global population. We therefore divide the "global" population into subpopulations representing sampled wastewater sites' catchment populations, with an additional subpopulation to represent
individuals who do not contribute to the sampled wastewater.

We model infection dynamics in these subpopulations hierarchically: subpopulation infection dynamics are distributed about a central jurisdiction-level infection dynamic, and the total infections that generate the hospital admissions observations are simply the sum of the subpopulation-level infections.

#### Subpopulation definition
The total population consists of $K_\mathrm{total}$ subpopulations $k$ with corresponding population sizes $n_k$. We associate one subpopulation to each of the $K_\mathrm{sites}$ wastewater sampling sites in the jurisdiction and assign that subpopulation a population size $n_k$ equal to the wastewater catchment population size reported for that wastewater catchment area.

Whenever the sum of the wastewater catchment population sizes $\sum\nolimits_{k=1}^{K_\mathrm{sites}} n_k$ is less than the total population size $n$, we use an additional subpopulation of size $n - \sum\nolimits_{k=1}^{K_\mathrm{sites}} n_k$ to model individuals in the population who are not covered by wastewater sampling.

The total number of subpopulations is then $K_\mathrm{total} = K_\mathrm{sites} + 1$: the $K_\mathrm{sites}$ subpopulations with sampled wastewater, and the final subpopulation to account for individuals not covered by wastewater sampling.

This amounts to modeling the wastewater catchments populations as approximately non-overlapping; every infected individual either does not contribute to measured wastewater or contributes principally to one wastewater catchment.
This approximation is reasonable if we restrict our analyses to primary wastewaster treatment plants, which avoids the possibility that an individual might be sampled once in a sample taken upstream and then sampled again in a more aggregated sample taken further downstream.

If the sum of the wastewater site catchment populations meets or exceeds the reported jurisdiction population ($\sum\nolimits_{k=1}^{K_\mathrm{sites}} n_k \ge n$) we do not use a final subpopulation without sampled wastewater. In that case, the total number of subpopulations $K_\mathrm{total} = K_\mathrm{sites}$.

When converting from predicted per capita incident hospital admissions $H(t)$ to predicted hospitalization counts, we use the jurisdiction population size $n$, even in the case where $\sum n_k > n$.

This amounts to making two key additional modeling assumptions:
- Any individuals who contribute to wastewaster measurements but are not part of the total population are distributed among the catchment populations approximately proportional to catchment population size.
- Whenever $\sum n_k \ge n$, the fraction of individuals in the jurisdiction not covered by wastewater is small enough to have minimal impact on the jurisdiction-wide per capita infection dynamics.

#### Subpopulation-level infections
We couple the subpopulation and total population infection dynamics at the level of the un-damped instantaneous reproduction number $\mathcal{R}^\mathrm{u}(t)$.

We model the subpopulations as having infection dynamics that are _similar_ to one another but can differ from the overall "global" dynamic.

We represent this with a hierarchical model where we first model a "global" un-damped effective reproductive number $\mathcal{R}^\mathrm{u}(t)$, but then allow individual subpopulations $k$ to have individual subpopulation values of $\mathcal{R}^\mathrm{u}_{k}(t)$

The "global" model for the undamped instantaneous reproductive number $\mathcal{R}^\mathrm{u}(t)$ follows the time-evolution described above.
Subpopulation deviations from the "global" reproduction number are modeled via a log-scale AR(1) process. Specifically, for subpopulation $k$:

$$
\log[\mathcal{R}^\mathrm{u}_{k}(t)] = \log[\mathcal{R}^\mathrm{u}(t)] + \delta_k(t)
$$

where $\delta_k(t)$ is the time-varying subpopulation effect on $\mathcal{R}(t)$, modeled as,

$$\delta_k(t) = \varphi_{R(t)} \delta_k(t-1) + \epsilon_{kt}$$

where $0 < \varphi_{R(t)} < 1$ and $\epsilon_{kt} \sim \mathrm{Normal}(0, \sigma_{R(t)\delta})$.

We chose a prior of $\varphi_{R(t)} \sim \mathrm{beta}(2,40)$ to impose limited autocorrelation in the week-by-week deviations.
We set a weakly informative prior $\sigma_{R(t)\delta} \sim \mathrm{Normal}(0, 0.3)$ to allow for either limited or substantial site-site heterogeneity in $\mathcal{R}(t)$, with the degree of heterogeneity inferred from the data.

The subpopulation $\mathcal{R}_{k}(t)$ is subject to the infection feedback described above such that:

```math
\mathcal{R}_k(t) = \mathcal{R}^\mathrm{u}_k(t) \exp \left(-\gamma \sum_{\tau = 1}^{T_f} I_k(t-\tau) g(\tau) \right)
```

From $\mathcal{R}_{k}(t)$, we generate values of the supopulation-level _expected_ latent incident infections per capita $I_k(t)$ using the renewal process described in [the infection component](#infection-component).

To obtain the number of infections per capita $I(t)$ in the total population as a whole, we sum the $K_\mathrm{total}$ subpopulation per capita infection counts $I_k(t)$ weighted by their population sizes:

```math
I(t) = \frac{1}{\sum\nolimits_{k=1}^{K_\mathrm{total}} n_k} \sum_{k=1}^{K_\mathrm{total}} n_k I_k(t)
```

We infer the site level initial per capita incidence $I_k(0)$ hierarchically. Specifically, we treat $\mathrm{logit}[I_k(0)]$ as Normally distributed about the corresponding jurisdiction-level value $\mathrm{logit}[I(0)]$, with an estimated standard deviation $\sigma_{i0}$:

$$
\mathrm{logit}[I_k(0)] \sim \mathrm{Normal}(\mathrm{logit}[I(0)], \sigma_{k0})
$$



### Hospital admissions component

Following other semi-mechanistic renewal frameworks, we model the _expected_ hospital admissions per capita $H(t)$ as a convolution of the _expected_ latent incident infections per capita $I(t)$, and a discrete infection to hospitalization distribution $d(\tau)$, scaled by the probability of being hospitalized $p_\mathrm{hosp}(t)$.

To account for day-of-week effects in hospital reporting, we use an estimated _weekday effect_ $\omega(t)$.
If $t$ and $t'$ are the same day of the week, $\omega(t) = \omega(t')$.
The seven values that $\omega(t)$ takes on are constrained to be non-negative and have a mean of 1.
This allows us to model the possibility that certain days of the week could have systematically high or low admissions reporting while holding the predicted weekly total reported admissions constant (i.e. the predicted weekly total is the same with and without these day-of-week reporting effects).

$$H(t) = \omega(t) p_\mathrm{hosp}(t) \sum_{\tau = 0}^{T_d} d(\tau) I(t-\tau)$$

Where $T_d$ is the maximum delay from infection to hospitalization that we consider.

We define the discrete hospital admissions delay distribution $d(\tau)$ as a convolution of the incubation period distribution [^Park2023] and a separate estimate of the distribution of time from symptom onset to hospital admission (see [Parameter section](#example-model-parameters-for-covid-19) below for further details).
#### Infection-hospitalization rate
In the models that include fits to wastewater data, we allow the population-level infection-hospitalization rate (IHR) to change over time. An inferred change in the IHR could reflect either a true change in the rate at which infections result in hospital admissions (e.g. the age distribution of cases could shift, a more or less severe variant could emerge, or vaccine coverage could shift) or a change in the relationship between infections and genomes shed in wastewater $G$ (which we currently treat as fixed, but which could change in time if, for example, immunity reduces wastewater shedding without reducing transmission, or a variant emerges with a different per infection wastewater shedding profile).

Therefore, we model the proportion of infections that give rise to hospital admissions $p_\mathrm{hosp}(t)$ as a piecewise-constant function with weekly change points.
If $t$ and $t'$ are two days in the same week, then $p_\mathrm{hosp}(t) = p_\mathrm{hosp}(t')$.

The values $p_\mathrm{hosp}(t)$ follow a logit-scale AR(1) process with linear scale median $\mu_{p_\mathrm{hosp}}$


$$ \mathrm{logit} (p_{\mathrm{hosp}}(t)) =  \mathrm{logit} (\mu_{p_{\mathrm{hosp}}}(t)) + \delta_{H}(t)$$

where $\delta_{H}(t)$ is the time-varying site effect on $\mathrm{logit} (p_{\mathrm{hosp}}(t))$, modeled as,

$$\delta_{H}(t) = \varphi_H \delta_{H}(t-1) + \epsilon_{Ht}$$

where $0 < \varphi_H < 1$ and $\epsilon_{Ht} \sim \mathrm{Normal}(0, \sigma_{H\delta})$.

We chose a relatively strong prior of $\varphi_H \sim \mathrm{beta}(1,100)$ to impose minimal autocorrelation, and similarly chose a small prior on $\sigma_{H\delta} \sim \mathrm{Normal}(0, 0.01)$ to impose our believe that the IHR should only change in time if the data strongly suggests it.
See [Prior Distributions](#prior-distributions) for the specified prior on $\mu_{p_\mathrm{hosp}}$.

In the version of the model where we do not fit the wastewater data (which we refer to as the "hospital admissiosn only" model), we model the IHR as constant. We assign this constant IHR the same prior distribution that we assign $\mu_{p_\mathrm{hosp}}$ in the wastewater-informed model.

#### Hospital admissions observation model
We model the observed hospital admission counts $h_t$ as:

$$h_t \sim \mathrm{NegBinom}(n H(t), \phi)$$

where the "global" population size $n$ (e.g. that of the state or county representing the catchment area population producing the hospital admissions) is used to convert from per-capita hospitalization rate $H(t)$ to hospitalization counts.

Currently, we do not explicitly model the delay from hospital admission to reporting of hospital admissions.
In reality, corrections (upwards or downwards) in the admissions data after the report date are possible and do happen.
This is an active area of further development, but for now, we advise the user to manually exclude hospital admissions data points that appear implausible.
Future work will include incorporation of a simple model for right-truncation when data is rolling in in real-time with incomplete reporting in recent days.
However, the current workflow assumes mandatory and for the most part complete reporting of hospital admissions.

### Viral genome concentration in wastewater component

We model site-specific viral genome concentrations in wastewater $C_i(t)$ independently for each site $i$, with the latent incident infections in subpopulation $k$ being mapped to the corresponding site $i$.
We model viral genome concentrations in wastewater in site $i$, $C_i(t)$ as a convolution of the _expected_ latent incident infections per capita in the corresponding subpopulation $k$, $I_k(t)$ and a normalized shedding kinetics function $s(\tau)$, multiplied by  $G$ the number of genomes shed per infected individual over the course of their infection and divided by $\alpha$ the volume of wastewater generated per person per day:

$$C_i(t) = \frac{G}{\alpha} \sum_{\tau = 0}^{\tau_\mathrm{shed}} s(\tau) I_k(t-\tau)$$

where $\tau_\mathrm{shed}$ is the total duration of fecal shedding and $i = k$ for all  $K_\mathrm{sites}$  wastewater sites.
Note there is no need to scale by wastewater catchment population size because $I_k(t)$ is measured as new infections per capita.

This approach assumes that $G$ and $\alpha$ are constant through time and across individuals.
In fact, there is substantial inter-individual variability in shedding kinetics and total shedding.
This approximation is more accurate when population sizes are large.
Incorporating the expected variability in the observed concentrations as a function of the number of contributing infected individduals is an area of active development.

We model the shedding kinetics $s(\tau)$ as a discretized, scaled triangular distribution[^Larremore2021]:

```math
\log_{10}[s^\mathrm{cont}(\tau)] = \begin{cases}
  V_\mathrm{peak} \frac{\tau}{\tau_\mathrm{peak}} & \tau \leq \tau_\mathrm{peak} \\
  V_\mathrm{peak} \left( 1 - \frac{\tau - \tau_\mathrm{peak}}{\tau_\mathrm{shed} - \tau_\mathrm{peak}} \right) & \tau_\mathrm{peak} < \tau \leq \tau_\mathrm{shed} \\
  0 & \tau > \tau_\mathrm{shed}
\end{cases}
```

where $V_\mathrm{peak}$ is the peak number or viral genomes shed on any day, $\tau_\mathrm{peak}$ is the time from infection to peak shedding, and $\tau_\mathrm{shed}$ is the total duration of shedding. Then:

Future iterations of this model will evaluate the utility of mechanistic modeling of wastewater collection and processing.

#### Viral genome concentration observation model

Genome concentration measurements can vary between sites, and even within a site through time, because of differences in sample collection and lab processing methods. To account for this variability, we add a scaling term $M_{ij}$ and a variablity term $\sigma_{cij}$ that vary across sites $i$ and also within sites across labs $j$:

$$\log[c_{ijt}] \sim \mathrm{Normal}(\log[M_{ij} C_i(t)], \sigma_{cij})$$

Both $M_{ij}$ and $\sigma_{cij}$ are modeled as site-level random effects.

In the rare cases when a site submits multiple concentrations for a single date and lab method, we treat each record as an independent observation.

### Censoring of wastewater observations below the limit of detection

Lab processing methods have a finite limit of detection (LOD), such that not all wastewater measurements can be modeled using the log-normal approach above.
This limit of detection varies across sites, between methods, and potentially also over time.

If an observed value $c_{ijt}$ is above the corresponding LOD, then the likelihood is:

$$
f_\mathrm{Normal}(\log[c_{ijt}]; \log[M_{ij} C_i(t)], \sigma_{cij})
$$

where $f_\mathrm{Normal}(x; \mu, \sigma)$ is the probability density function of the Normal distribution.
When the observed value is below the LOD, we use a censored likelihood:

```math
\int_{-\infty}^{\log [\mathrm{LOD}_{ijt}]} f_\mathrm{Normal}(x; \log[M_{ij} C_i(t)], \sigma_{cij}) \mathrm{d}x
```

(This is mathematically equivalent to integrating the probability density function of the log-normal distribution from zero to the LOD.)


## Example model parameters for COVID-19

The default parameters provided by the `wwinference` package are used to fit a model of COVID-19 hospital admissions and wastewater concentrations in terms of reported SARS-CoV-2 genome copies per mL.
Below we will describe the priors and parameters provided. If fitting the model to a different epidemiological indicator (e.g. cases) or a different pathogen (e.g. flu) a number of these will have to be modified accordingly.

### Prior distributions

We use informative priors for parameters that have been well characterized in the literature and weakly informative priors for parameters that have been less well characterized.

| Parameter | Prior distribution | Source |
|---|---|---|
| Initial hospitalization probability | $\mathrm{logit}[p_{\mathrm{hosp}}(t_0)] \sim \mathrm{Normal}(\mathrm{logit}[0.01], 0.3)$ | Perez-Guzman et al. 2023 [^Perez] |
| Time to peak fecal shedding |  $\tau_\mathrm{peak} \sim \mathrm{Normal}(5 \text{ days}, 1 \text{ day})$ | Russell et al. 2023 [^Russell], Huisman et al. 2022 [^Huisman], Cavany et al. 2022 [^Cavany] |
| Peak viral shedding $V_\mathrm{peak}$| $\log_{10}[V_\mathrm{peak}] \sim \mathrm{Normal}(5.1, 0.5)$ | Miura et al. 2021 [^Muira] |
| Duration of shedding | $\tau_\mathrm{shed} \sim \mathrm{Normal}(17 \text{ days}, 3 \text{ days})$  | Cevik et al. 2021 [^Cevik], Russell et al. 2023 [^Russell]   |
| Total genomes shed per infected individual | $\log_{10}[G] \sim \mathrm{Normal}(9, 2)$ | Watson et al 2023[^Watson]   |
| Initial infections per capita $I_0$ | $I_0 \sim \mathrm{Beta}(1 + k i_\mathrm{est}, 1 + k (1-i_\mathrm{est}))$ | where $i_\mathrm{est}$ is the sum of the last 7 days of hospital admissions, divided by jurisdiction population, and divided by the prior mode for $p_\mathrm{hosp}$, and $k = 5$ is a parameter governing the informativeness ("certainty") of the Beta distribution |
| Initial exponential growth rate | $r \sim \mathrm{Normal}(0, 0.01)$ | Chosen to assume flat dynamics prior to observations |
| Infection feedback term | $\gamma \sim \mathrm{logNormal}(6.37, 0.4)$ | Weakly informative prior chosen to have a mode of 500 in natural scale, based on posterior estimates of peaks from prior seasons in a few jurisdictions |

### Scalar parameters

| Parameter | Value | Source |
|---|---|---|
| Maximum generation interval | $T_g = 15$ days | |
| Maximum infection to hospital admissions delay | $T_d = 55$ days| |
| Wastewater produced per person-day | $\alpha=$ 378,500 mL per person-day | Ortiz 2024[^Ortiz] |

### Distributional parameters

The discrete generation interval probability mass function $g(\tau)$ approximates a log-normal distribution[^Park2023] with log-mean 2.9 and log-standard deviation of 1.64.
To approximate the double censoring process necessary to discretize the continuous log-normal distribution, we use a simulation-based approach as recommended by Park et al.[^Park2024].
This assumes that the primary event is uniformly distributed (this ignores the influence of the growth rate within the primary interval but is a good approximation in most settings). The secondary event is then a sum of this primary interval and the continuous distribution and is observed within a day (see Figure 9 in [^Park2024]).
As the renewal process is not defined if there is probability mass on day zero we further left truncate this distribution.
For more details refer to [^Park2024].

We derive the distribution $\delta(\tau)$ of the delay from infection to hospital admission as the sum of the incubation period (delay from infection to symptom onset) and the period from symptom onset to hospital admission.

We model the incubation period with a discretized, modified Weibull distribution[^Park2023] with probability mass function $\delta(\tau)$:

$$
\delta^\mathrm{cont}(\tau) = \exp[0.15\tau] f_\mathrm{Weibull}(\tau; \mathrm{shape}=1.5, \mathrm{scale}=3.6)
$$

```math
\delta(\tau) = \begin{cases}
\delta^\mathrm{cont}(\tau) / \left( {\sum}_{\tau'=0}^{23} g^\mathrm{cont}(\delta') \right) & 0 \leq \tau \leq 23 \\
0 & \text{otherwise}
\end{cases}
```

We model the symptom onset to hospital admission delay distribution with a Negative Binomial distribution with probability mass function $\gamma(\tau)$ fit to line list patient data from Dananché et al. 2022[^Danache2022].

$$
\gamma(\tau) = f_\mathrm{NegBin}(\tau; 6.99 \text{ days}, 2.49 \text{ days})
$$

The infection-to-hospitalization delay distribution $d(\tau)$ is the convolution:

$$
d(\tau) = \sum_{x=0}^\tau \delta(x) \gamma(\tau - x)
$$

This resulting infection to hospital admission delay distribution has a mean of 12.2 days and a standard deviation of 5.67 days.

## Implementation

Our framework is an extension of the widely used [^CDCRtestimates] [^CDCtechnicalblog], semi-mechanistic renewal framework `{EpiNow2}` [^epinow2paper][^EpiNow2], using a Bayesian latent variable approach implemented in the probabilistic programming language Stan [^stan] using [^cmdstanr] to interface with R.

We fit the model using the “No-U-Turn Sampler Markov chain Monte Carlo” method. This is a type of Hamiltonian Monte Carlo (HMC) algorithm and is the core fitting method used by `cmdstan`.
The dfault parameter settings are set to run 4 chains for 750 warm-up iterations and 500 sampling iterations, with a target average acceptance probability of 0.95 and a maximum tree depth of 12.
The user can adjust these settings using the `get_mcmc_options()` function.

## Appendix: Wastewater data pre-processing

### Viral genome concentration in wastewater outlier detection and removal

We identify potential outlier genome concentrations for each unique site and lab pair with an approach based on $z$-scores.

Briefly, we compute $z$-scores for the concentrations and their finite differences and remove any observations above a threshold values for either metric. In detail:

1. For purposes of outlier detection, exclude wastewater observations below the LOD.
1. For purposes of outlier detection, exclude observations more than 90 days before the forecast date.
1. For each site $i$, compute the change per unit time between successive observations $t$ and $t'$: $(\log[c_{it'}] - \log[c_{it}])/(t' - t)$.
1. Compute $z$-scores for $\log[c_{it}]$ across all sites $i$ and timepoints $t$. Flag values with $z$-scores over 3 as outliers and remove them from model calibration.
1. Compute $z$-scores for the change per unit time values across all sites and pairs of timepoints. For values with $z$-scores over 2, flag the corresponding wastewater concentrations $c_{it}$ as outliers and remove them from model calibration.

The $z$-score thresholds were chosen by visual inspection of the data.

## Appendix: Notation

The notation $X \sim \mathrm{Distribution}$ indicates that a random variable $X$ is distributed according to a given distribution.

We parameterize Normal distributions in terms of their mean and standard deviation: $\mathrm{Normal}(\mathrm{mean, standard\ deviation})$.

We parameterize Beta distributions in terms of their two standard shape parameters $\alpha$ and $\beta$, which can be [interpreted in terms of the counts of observed successes and failures](https://stats.stackexchange.com/questions/47771/what-is-the-intuition-behind-beta-distribution), respectively, in a Binomial experiment to infer a probability: $\mathrm{Beta}(\alpha, \beta)$.

We parameterize Negative Binomial distributions in terms of their mean and their positive-constrained dispersion parameter (often denoted $\phi$): $\mathrm{NegBinom}(\mathrm{mean, dispersion})$. As the dispersion parameter goes to 0, a Negative Binomial distribution becomes increasingly over-dispersed. As it goes to positive infinity, the Negative Binomial approximates a Poisson distribution with the same mean.

We write $\mathrm{logit}(x)$ to refer to the logistic transform: $\mathrm{logit}(x) \equiv \log(x) - \log(1 - x)$.

Observed data are labeled by data source: $c$ for wastewater concentrations, $h$ for hospital admissions.
Hospitalization data are indexed by day $t$ (i.e., $h_t$).
Wastewater data are indexed by site $i$, wastewater testing lab $j$, and day $t$ (e.g., $c_{ijt}$).


## References

[^epinow2paper]: Abbott, S. et al. Estimating the time-varying reproduction number of SARS-CoV-2 using national and subnational case counts. _Wellcome Open Res_. 5:112 (2020). https://doi.org/10.12688/wellcomeopenres.16006.2
[^Asher2018]: Asher, J. Forecasting Ebola with a regression transmission model. _Epidemics._ **22**, 50-55 (2018). https://doi.org/10.1016/j.epidem.2017.02.009
[^stan]: Stan Development Team. _Stan Modeling Language Users Guide and Reference Manual_. (2023). https://mc-stan.org
[^CDCRtestimates]: US Centers for Disease Control and Prevention. _Current Epidemic Growth Status (Based on Rt) for States and Territories_. https://www.cdc.gov/forecast-outbreak-analytics/about/rt-estimates.html (2024).
[^CDCtechnicalblog]: US Centers for Disease Control and Prevention. _Technical Blog: Improving CDC’s Tools for Assessing Epidemic Growth_ (2024). https://www.cdc.gov/forecast-outbreak-analytics/about/technical-blog-rt.html
[^Park2023]: Park, S.W. et al. Inferring the differences in incubation-period and generation-interval distributions of the Delta and Omicron variants of SARS-CoV-2. _Proc Natl Acad Sci U S A_. 120(22):e2221887120 (2023). https://doi.org/10.1073/pnas.2221887120
[^Danache2022]: Danaché, C. et al. Baseline clinical features of COVID-19 patients, delay of hospital admission and clinical outcome: A complex relationship. _PLoS One_ 17(1):e0261428 (2022). https://doi.org/10.1371/journal.pone.0261428
[^Perez]: Perez-Guzman, P.N. et al. Epidemiological drivers of transmissibility and severity of SARS-CoV-2 in England. _Nat Commun_ 14, 4279 (2023). https://doi.org/10.1038/s41467-023-39661-5
[^Muira]: Miura F, Kitajima M, Omori R. Duration of SARS-CoV-2 viral shedding in faeces as a parameter for wastewater-based epidemiology: Re-analysis of patient data using a shedding dynamics model. _Sci Total Environ_ 769:144549 (2021). https://doi.org/10.1016/j.scitotenv.2020.144549
[^Huisman]: Huisman, J.S. et al. Estimation and worldwide monitoring of the effective reproductive number of SARS-CoV-2 _eLife_ 11:e71345 (2022). https://doi.org/10.7554/eLife.71345
[^Cavany]: Cavany S, et al. Inferring SARS-CoV-2 RNA shedding into wastewater relative to the time of infection. _Epidemiology and Infection_ 150:e21 (2022). https://doi.org/10.1017/S0950268821002752
[^Russell]: Russell, T.W. et al. Within-host SARS-CoV-2 viral kinetics informed by complex life course exposures reveals different intrinsic properties of Omicron and Delta variants. _medRxiv_ (2023).  https://doi.org/10.1101/2023.05.17.23290105
[^Cevik]: Cevik, M. et al. SARS-CoV-2, SARS-CoV, and MERS-CoV viral load dynamics, duration of viral shedding, and infectiousness: a systematic review and meta-analysis. _Lancet Microbe_ **2(1)**,e13-e22 (2021). https://doi.org/10.1016/S2666-5247(20)30172-5
[^Watson]: Leighton, M. et al. Improving estimates of epidemiological quantities by combining reported cases with wastewater data: a statistical framework with applications to COVID-19 in Aotearoa New Zealand. _medRxiv_ (2023). https://doi.org/10.1101/2023.08.14.23294060
[^Ortiz]: Ortiz, P. _Wastewater facts - statistics and household data in 2024_. https://housegrail.com/wastewater-facts-statistics/
[^Larremore2021]: Larremore, D.B. et al. Test sensitivity is secondary to frequency and turnaround time for COVID-19 screening. _Science Advances_ (2021). https://doi.org/10.1126/sciadv.abd5393
[^Cori]: Cori, A., Ferguson, N. M., Fraser, C., & Cauchemez, S. A new framework and software to estimate time-varying reproduction numbers during epidemics. _Am. J. Epidemiol._ **178**, 1505-1512 (2013). https://doi.org/10.1093/aje/kwt133
[^EpiNow2]: Abbott, S. et al. _EpiNow2: Estimate real-time case counts and time-varying epidemiological parameters._ https://doi.org/10.5281/zenodo.3957489
[^Epidemia]: Fraser, C. (2007). Estimating individual and household reproduction numbers in an emerging epidemic. _PLoS One_, **2**(8), e758 (2007). https://doi.org/10.1371/journal.pone.0000758
[^cmdstanr]: _CmdStanR: the R interface to CmdStan_. (2024). https://mc-stan.org/cmdstanr/index.html
[^Park2024]: Park, S.W. et al. Estimating epidemiological delay distributions for infectious diseases.
_medRxiv_ (2024). https://doi.org/10.1101/2024.01.12.24301247
[^Gostic2020]: Gostic, K.M. et al. Practical Considerations for Measuring the Effective Reproductive Number, Rt. _PLoS Comput Biol_. **16**(12) (2020). https://doi.org/10.1371/journal.pcbi.1008409
