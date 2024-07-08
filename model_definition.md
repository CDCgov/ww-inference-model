## Model overview: wastewater-informed, site-level infection dynamics 

This document details the generative model used by `wwinference` to infer global and local infection dynamics from count data (e.g. cases or hospital admissions) and wastewater concentration data and use that to produce nowcasts and forecasts. 
The `wwinference` model assumes that wastewater concentration data is available for one or more "local" sites that represent a subset of the total population that produce the "global" epidemiological indicator (cases or admissions).
In the future, we plan to provide the user with functionality for other types of data structures, e.g. multiple data streams of hospital admissions and multiple wastewater concentration data streams.
The model is structured into three models: infections(#infection-component), hospital admissions(#hospital-admissions-component), and viral genome concentrations in wastewater(#wastewater-component)

### Model components

Our models are constructed from a set of generative components. These are:

- [**Infection component:**](#infection-component) A renewal model for the infection dynamics, which generates estimates of incident latent infections per capita.
- [**Hospital admissions component:**](#hospital-admissions-component) A model for the expected number of hospital admissions given incident latent infections per capita.
- [**Viral genome concentration in wastewater:**](#wastewater-component) A model for the expected genome concentration given incident infections per capita.

Depending on the model, these components are implemented at different spatial scales and with different observation processes. In particular, the link to the observables depends on the model and the form of the observables, see below for detailed descriptions of each component and the following sections for model-specific details.

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
\log[\mathcal{R}^\mathrm{u}(t_3)] \sim \mathrm{Normal}\left(\log[\mathcal{R}^\mathrm{u}(t_2)] + \beta \left(\log[\mathcal{R}^\mathrm{u}(t_2)] - \log[\mathcal{R}^\mathrm{u}(t_1)]\right), \sigma_r \right)
$$

where $t_1$, $t_2$, and $t_3$ are days in three successive weeks, $\beta$ is an autoregression coefficient which serves to make week-to-week changes correlated, and $\sigma_r$ determines the overall variation in week-to-week changes.
We bound $\beta$ to be between 0 and 1 so that any changes in trend in $\mathcal{R}^\mathrm{u}(t)$ are damped over time to a degree determined by $\beta$.

The damping term we use is based on Asher et al. 2018[^Asher2018] but extended to be applicable to a renewal process. It assumes that the instantaneous reproduction number is damped by recent infections weighted by the generation interval. This is a simple way to account for the fact that the instantaneous reproduction number is likely to decrease when there are many infections in the population, due to factors such as immunity, behavioral changes, and public health interventions. The damping term is defined as:

$$ \mathcal{R}(t) = \mathcal{R}^\mathrm{u}(t) \exp \left( -\gamma \sum_{\tau = 1}^{T_f}I(t-\tau)g(\tau) \right) $$

where $\gamma$ is the _infection feedback term_ controlling the strength of the damping on $\mathcal{R}(t)$, and the summation is analogous to the "force of infection." See [Prior Distributions](#prior-distributions) below for description of prior choice on $\gamma$.

### Hospital admissions component

Following other semi-mechanistic renewal frameworks, we model the _expected_ hospital admissions per capita $H(t)$ as a convolution of the _expected_ latent incident infections per capita $I(t)$, and a discrete infection to hospitalization distribution $d(\tau)$, scaled by the probability of being hospitalized $p_\mathrm{hosp}(t)$.

To account for day-of-week effects in hospital reporting, we use an estimated _weekday effect_ $\omega(t)$. If $t$ and $t'$ are the same day of the week, $\omega(t) = \omega(t')$.
The seven values that $\omega(t)$ takes on are constrained to be non-negative and have a mean of 1.
This allows us to model the possibility that certain days of the week could have systematically high or low admissions reporting while holding the predicted weekly total reported admissions constant (i.e. the predicted weekly total is the same with and without these day-of-week reporting effects).

$$H(t) = \omega(t) p_\mathrm{hosp}(t) \sum_{\tau = 0}^{T_d} d(\tau) I(t-\tau)$$

Where $T_d$ is the maximum delay from infection to hospitalization that we consider.

We define the discrete hospital admissions delay distribution $d(\tau)$ as a convolution of the incubation period distribution [^Park2023] and a separate estimate of the distribution of time from symptom onset to hospital admission (see [Parameter section](#model-parameters) below for further details).
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

In hospital admissions only models, we model the IHR as constant. We assign this constant IHR the same prior distribution that we assign $\mu_{p_\mathrm{hosp}}$ in the wastewater model.

### Hospital admissions observation model
We model the observed hospital admission counts $h_t$ as:

$$h_t \sim \mathrm{NegBinom}(n H(t), \phi)$$

where the jurisdiction population size $n$ is used to convert from per-capita hospitalization rate $H(t)$ to hospitalization counts.

Currently, we do not explicitly model the delay from hospital admission to reporting of hospital admissions. In reality, corrections (upwards or downwards) in the admissions data after the report date are possible and do happen. See [outlier detection and removal](#appendix-wastewater-data-pre-processing) for further details.

### Viral genome concentration in wastewater component

We model viral genome concentrations in wastewater $C(t)$ as a convolution of the _expected_ latent incident infections per capita $I(t)$ and a normalized shedding kinetics function $s(\tau)$, multiplied by  $G$ the number of genomes shed per infected individual over the course of their infection and divided by $\alpha$ the volume of wastewater generated per person per day:

$$C(t) = \frac{G}{\alpha} \sum_{\tau = 0}^{\tau_\mathrm{shed}} s(\tau) I(t-\tau)$$

where $\tau_\mathrm{shed}$ is the total duration of fecal shedding.
Note there is no need to scale by wastewater catchment population size because $I(t)$ is measured as new infections per capita.

This approach assumes that $G$ and $\alpha$ are constant through time and across individuals.
In fact, there is substantial inter-individual variability in shedding kinetics and total shedding.
This approximation is more accurate when population sizes are large.

We model the shedding kinetics $s(\tau)$ as a discretized, scaled triangular distribution[^Larremore2021]:

```math
\log_{10}[s^\mathrm{cont}(\tau)] = \begin{cases}
  V_\mathrm{peak} \frac{\tau}{\tau_\mathrm{peak}} & \tau \leq \tau_\mathrm{peak} \\
  V_\mathrm{peak} \left( 1 - \frac{\tau - \tau_\mathrm{peak}}{\tau_\mathrm{shed} - \tau_\mathrm{peak}} \right) & \tau_\mathrm{peak} < \tau \leq \tau_\mathrm{shed} \\
  0 & \tau > \tau_\mathrm{shed}
\end{cases}
```

where $V_\mathrm{peak}$ is the peak number or viral genomes shed on any day, $\tau_\mathrm{peak}$ is the time from infection to peak shedding, and $\tau_\mathrm{shed}$ is the total duration of shedding. Then:

We model the log observed genome concentrations as Normally distributed:

$$
\log[c_t] \sim \mathrm{Normal}(C(t), \sigma_c)
$$

This component does not mechanistically simulate each step involved in sample collection, processing, and reporting.
Instead, it aims to to account for these processes, at a summary level.
Future iterations of this model will evaluate the utility of mechanistic modeling of wastewater collection and processing.

## Model 1: Site-level infection dynamics

In this model, we represent hospital admissions at the jurisdictional level and viral genome concentrations at the site level. We use the components described above but divide the jurisdiction's total population into subpopulations representing sampled wastewater sites' catchment populations, with an additional subpopulation to represent individuals who do not contribute to sampled wastewater.

We model infection dynamics in these subpopulations hierarchically: subpopulation infection dynamics are distributed about a central jurisdiction-level infection dynamic, and jurisdiction's total infections are simply the sum of the subpopulation-level infections.

### Subpopulation definition
In Model 1, a jurisdiction consists of $K_\mathrm{total}$ subpopulations $k$ with corresponding population sizes $n_k$. We associate one subpopulation to each of the $K_\mathrm{sites}$ wastewater sampling sites in the jurisdiction and assign that subpopulation a population size $n_k$ equal to the wastewater catchment population size reported to NWSS for that site.

Whenever the sum of the wastewater catchment population sizes $\sum\nolimits_{k=1}^{K_\mathrm{sites}} n_k$ is less than the total jurisdiction population size $n$, we use an additional subpopulation of size $n - \sum\nolimits_{k=1}^{K_\mathrm{sites}} n_k$ to model individuals in the jurisdiction who are not covered by wastewater sampling.

The total number of subpopulations is then $K_\mathrm{total} = K_\mathrm{sites} + 1$: the $K_\mathrm{sites}$ subpopulations with sampled wastewater, and the final subpopulation to account for individuals not covered by wastewater sampling.

This amounts to modeling the wastewater catchments populations as approximately non-overlapping; every infected individual either does not contribute to measured wastewater or contributes principally to one wastewater catchment.
This approximation is reasonable because we only use samples taken from primary wastewaster treatment plants, which avoids the possibility that an individual might be sampled once in a sample taken upstream and then sampled again in a more aggregated sample taken further downstream; see [data filtering](#appendix-wastewater-data-pre-processing) for further details.

If the sum of the wastewater site catchment populations meets or exceeds the reported jurisdiction population ($\sum\nolimits_{k=1}^{K_\mathrm{sites}} n_k \ge n$) we do not use a final subpopulation without sampled wastewater. In that case, the total number of subpopulations $K_\mathrm{total} = K_\mathrm{sites}$. In our data, this is true only for the District of Columbia.

When converting from predicted per capita incident hospital admissions $H(t)$ to predicted hospitalization counts, we use the jurisdiction population size $n$, even in the case of the District of Columbia where $\sum n_k > n$.

This amounts to making two key additional modeling assumptions:
- Any individuals who contribute to wastewaster measurements but are not part of the jurisdiction population are distributed among the catchment populations approximately proportional to catchment population size.
- Whenever $\sum n_k \ge n$, the fraction of individuals in the jurisdiction not covered by wastewater is small enough to have minimal impact on the jurisdiction-wide per capita infection dynamics.

### Subpopulation-level infection dynamics
We couple the subpopulation- and jurisdiction-level infection dynamics at the level of the un-damped instantaneous reproduction number $\mathcal{R}^\mathrm{u}(t)$.

We model the subpopulations as having infection dynamics that are _similar_ to one another but can differ from the overall jurisdiction-level dynamic.

We represent this with a hierarchical model where we first model a jurisdiction-level un-damped effective reproductive number $\mathcal{R}^\mathrm{u}(t)$, but then allow individual subpopulations $k$ to have individual subpopulation values of $\mathcal{R}^\mathrm{u}_{k}(t)$

The jurisdiction-level model for the undamped instantaneous reproductive number $\mathcal{R}^\mathrm{u}(t)$ follows the time-evolution described above.
Subpopulation deviations from the jurisdiction-level reproduction number are modeled via a log-scale AR(1) process. Specifically, for subpopulation $k$:

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

To obtain the number of infections per capita $I(t)$ in the jurisdiction as a whole, we sum the $K_\mathrm{total}$ subpopulation per capita infection counts $I_k(t)$ weighted by their population sizes:

```math
I(t) = \frac{1}{\sum\nolimits_{k=1}^{K_\mathrm{total}} n_k} \sum_{k=1}^{K_\mathrm{total}} n_k I_k(t)
```

We infer the site level initial per capita incidence $I_k(0)$ hierarchically. Specifically, we treat $\mathrm{logit}[I_k(0)]$ as Normally distributed about the corresponding jurisdiction-level value $\mathrm{logit}[I(0)]$, with an estimated standard deviation $\sigma_{i0}$:

$$
\mathrm{logit}[I_k(0)] \sim \mathrm{Normal}(\mathrm{logit}[I(0)], \sigma_{k0})
$$

### Viral genome concentration in wastewater

We model site-specific viral genome concentrations in wastewater $C_i(t)$ independently for each site $i$ using the same model as described in [the wastewater component](#wastewater-viral-concentration-component). The latent incident infections in subpopulation $k$ are mapped to the corresponding site $i$.

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

If a sample is flagged in the NWSS data as below the LOD (field `pcr_target_below_lod`) but is missing a reported LOD (field `lod_sewage`), the 95th percentile of LOD values across the entire data is used as the integral's upper limit.

If a sample has a reported concentration (field `pcr_target_avg_conc`) above the corresponding reported LOD, but the sample is nevertheless flagged as below the LOD (field `pct_target_below_lod`), we assume the flag takes precedence and treat the sample as below LOD for the purposes of censoring.
