functions {
#include functions/ar1.stan
#include functions/diff_ar1.stan
#include functions/convolve.stan
#include functions/hospitalization.stan
#include functions/infections.stan
#include functions/observation_model.stan
#include functions/utils.stan
#include functions/construct_spatial_rt.stan
#include functions/dist_matrix_normalization.stan
#include functions/matrix_normalization.stan
#include functions/independence_corr_func.stan
#include functions/exponential_decay_corr_func.stan
#include functions/construct_aux_rt.stan

}

// The fixed input data
data {
  int<lower=1> gt_max;
  int<lower=1> hosp_delay_max;
  vector<lower=0,upper=1>[hosp_delay_max] inf_to_hosp; // delay distribution from infecion to hospital admission
  real<lower=0> mwpd; // mL of ww produced per person per day
  int<lower=1> if_l; // length of infection feedback pmf
  vector<lower=0,upper=1>[if_l] infection_feedback_pmf; // infection feedback pmf
  int<lower=0> ot; // maximum time index for the hospital admissions (max number of days we could have observations)
  int<lower=0> oht; // number of days that we have hospital admissions observations
  int<lower=0> n_subpops; // number of WW sites
  int<lower=0> n_ww_lab_sites; // number of unique ww-lab combos
  int<lower=0> n_censored; // numer of observed WW data points that are below the LOD
  int<lower=0> n_uncensored; //number not below LOD
  int<lower=0> owt; // number of days of observed WW (should be roughly ot/7)
  int<lower=0> uot; // unobserved time before we observe hospital admissions/ WW
  int<lower=0> ht; // horizon time (nowcast + forecast time)
  int<lower=0> n_weeks; // number of weeks for weekly random walk on R(t)
  matrix<lower=0> [ot+ht, n_weeks] ind_m; // matrix to convert R(t) from weekly to daily
  int<lower=0> tot_weeks; // number of weeks for the weekly random walk on IHR (includes unobserved time)
  matrix<lower=0> [uot+ot+ht, tot_weeks] p_hosp_m; // matrix to convert p_hosp from weekly to daily
  vector<lower=0,upper=1>[gt_max] generation_interval; // generation interval distribution
  real<lower = 1e-20> state_pop; // population size
  vector<lower = 1e-20>[n_subpops] subpop_size; // the population sizes for each subpopulation
  real<lower = state_pop> norm_pop;
  array[owt] int<lower=1,upper=ot + ht> ww_sampled_times; // a list of all of the days on which WW is sampled
                                   // will be mapped to the corresponding sites (ww_sampled_sites)
  array[oht] int<lower=1, upper=ot> hosp_times; // the days on which hospital admissions are observed
  array[owt] int<lower=1,upper=n_subpops> ww_sampled_sites; // vector of unique sites in order of the sampled times
  array[owt] int<lower=1,upper=n_ww_lab_sites> ww_sampled_lab_sites; // vector of unique lab-site combos i
   	     // n order of the sampled times
  array[n_censored] int<lower=1,upper=owt> ww_censored; // times that the WW data is below the LOD
  array[n_uncensored] int<lower=1,upper=owt> ww_uncensored; // time that WW data is above LOD
  vector[owt] ww_log_lod; // The limit of detection in that site at that time point
  array[n_ww_lab_sites] int<lower=1,upper=n_subpops> lab_site_to_site_map; // which lab sites correspond to which sites
  array[oht] int<lower=0> hosp; // observed hospital admissions
  array[ot + ht] int<lower=1,upper=7> day_of_week; // integer vector with 1-7 corresponding to the weekday
  vector[owt] log_conc; // observed concentration of viral genomes in WW
  int<lower=0,upper=1> compute_likelihood; // 1= use data to compute likelihood
  int<lower=0,upper=1> include_ww; // 1= include wastewater data in likelihood calculation
  int<lower=0,upper=1> include_hosp; // 1 = fit to hosp, 0 = only fit wastewater model
  vector[6] viral_shedding_pars;// tpeak, viral peak, shedding duration mean and sd
  real<lower=0> autoreg_rt_a;
  real<lower=0> autoreg_rt_b;
  real<lower=0> autoreg_rt_site_a;
  real<lower=0> autoreg_rt_site_b;
  real<lower=0> autoreg_p_hosp_a;
  real<lower=0> autoreg_p_hosp_b;
  real inv_sqrt_phi_prior_mean;
  real<lower=0> inv_sqrt_phi_prior_sd;
  real r_prior_mean;
  real<lower=0> r_prior_sd;
  real log10_g_prior_mean;
  real<lower=0> log10_g_prior_sd;
  real<lower=0> i_first_obs_over_n_prior_a;
  real<lower=0> i_first_obs_over_n_prior_b;
  real sigma_i_first_obs_prior_mode;
  real<lower=0> sigma_i_first_obs_prior_sd;
  vector<lower=0>[7] hosp_wday_effect_prior_alpha;
  real mean_initial_exp_growth_rate_prior_mean;
  real<lower=0> mean_initial_exp_growth_rate_prior_sd;
  real sigma_initial_exp_growth_rate_prior_mode;
  real<lower=0> sigma_initial_exp_growth_rate_prior_sd;
  real mode_sigma_ww_site_prior_mode;
  real<lower=0> mode_sigma_ww_site_prior_sd;
  real sd_log_sigma_ww_site_prior_mode;
  real<lower=0> sd_log_sigma_ww_site_prior_sd;
  real<lower=0> eta_sd_sd;
  real p_hosp_prior_mean;
  real<lower=0> p_hosp_sd_logit;
  real<lower=0> p_hosp_w_sd_sd;
  real<lower=0> ww_site_mod_sd_sd;
  real<lower=0> sigma_rt_prior;
  real log_phi_g_prior_mean;
  real<lower=0> log_phi_g_prior_sd;
  real inf_feedback_prior_logmean;
  real<lower=0> inf_feedback_prior_logsd;

  real log_phi_mu_prior;
  real log_phi_sd_prior;
  real l;
  real log_sigma_generalized_mu_prior;
  real log_sigma_generalized_sd_prior;
  real log_scaling_factor_mu_prior;
  real log_scaling_factor_sd_prior;
  matrix[n_subpops-1, n_subpops-1] dist_matrix;
  int<lower=0, upper=2> corr_structure_switch;
}

// The transformed data
transformed data {
  // viral shedding parameters
  real t_peak_mean = viral_shedding_pars[1];
  real<lower=0> t_peak_sd = viral_shedding_pars[2];
  real viral_peak_mean = viral_shedding_pars[3];
  real<lower=0> viral_peak_sd = viral_shedding_pars[4];
  real dur_shed_mean = viral_shedding_pars[5];
  real<lower=0> dur_shed_sd = viral_shedding_pars[6];

  // natural scale -> lognormal parameters
  // https://en.wikipedia.org/wiki/Log-normal_distribution
  real r_logmean = convert_to_logmean(r_prior_mean, r_prior_sd);
  real<lower=0> r_logsd = convert_to_logsd(r_prior_mean, r_prior_sd);
  // reversed generation interval
  vector<lower=0,upper=1>[gt_max] gt_rev_pmf = reverse(generation_interval);
  vector<lower=0,upper=1>[if_l] infection_feedback_rev_pmf = reverse(infection_feedback_pmf);

  // normalizing dist matrix using largest dist
  matrix[n_subpops - 1, n_subpops - 1] norm_dist_matrix = dist_matrix_normalization(dist_matrix);
}

// The parameters accepted by the model.
parameters {
  vector[n_weeks-1] w; // weekly random walk of state-level mean baseline R(t) (log scale)
  real<lower=0> eta_sd;
  real<lower=0, upper=1> autoreg_rt;// coefficient on AR process in R(t)
  real log_r_mu_intercept; // state-level mean baseline reproduction number estimate (log) at t=0
  real<lower=0> sigma_rt; // magnitude of site level variation from state level
  real<lower=0, upper=1> autoreg_rt_site;
  real<lower=0, upper=1> autoreg_p_hosp;
  real<lower=0,upper=1> i_first_obs_over_n; // per capita
  // infection incidence on the day of the first observed infection
  vector[n_subpops] eta_i_first_obs; // z-score on logit scale of site
  // initial per capita infection incidence relative to state value
  real<lower=0> sigma_i_first_obs; // stdev between logit state and site initial
  // per capita infection incidence
  vector[n_subpops] eta_initial_exp_growth_rate; // z scores of individual site level initial exponential growth rates
  real<lower=0> sigma_initial_exp_growth_rate; // sd of distribution of site level initial exp growth rates
  real mean_initial_exp_growth_rate; // mean of distribution of site level initial exp growth rates
  real<lower=1/sqrt(5000)> inv_sqrt_phi_h;
  real mode_sigma_ww_site; //mode of site level stdev
  real<lower=0> sd_log_sigma_ww_site; // stdev of the log site level stdev
  vector[n_ww_lab_sites] eta_log_sigma_ww_site; // let each lab-site combo have its own observation error
  real p_hosp_mean; // Estimated mean IHR
  vector[tot_weeks] p_hosp_w; // weekly random walk for IHR
  real<lower=0> p_hosp_w_sd; // Estimated IHR sd
  real<lower=0> t_peak; // time to viral load peak in shedding
  real viral_peak; // log10 peak viral load shed /mL
  real<lower=0> dur_shed; // duration of detectable viral shedding
  real log10_g; // mean log10 of number of genomes per infected individual
  vector[n_ww_lab_sites] ww_site_mod_raw; // lab-site specific WW modifier on the observation error
  real<lower=0> ww_site_mod_sd; // site specific WW modifier stdev
                                // for now assumes the same across sites, can change or throw into
                                // observation error
  simplex[7] hosp_wday_effect; // day of week reporting effect, sums to 1
  real<lower=0> infection_feedback; // infection feedback

  // Site spatial params--------------------------------------------------------
  //matrix[n_subpops, n_subpops] non_norm_omega;
  real log_sigma_generalized;
  real log_phi;
  real log_scaling_factor;
  matrix[n_subpops-1,n_weeks] non_cent_spatial_dev_ns_mat;
  vector[n_weeks] norm_vec_aux_site;
  cholesky_factor_corr[corr_structure_switch == 2 ? n_subpops-1 : 2] L_Omega;
  //----------------------------------------------------------------------------
}
//
transformed parameters {
  vector[ot + uot + ht] p_hosp; // probability of hospitalization
  vector[ot + uot + ht] model_hosp_per_capita; // model estimated hospital admissions per capita
  vector[oht] exp_obs_hosp; //  expected observed hospital admissions
  vector[ot] exp_obs_hosp_per_capita_no_wday_effect; // expected observed hospital admissions per capita before weekday effect
  vector[gt_max] s; // viral kinetics trajectory (normalized)
  vector[owt] exp_obs_log_v_true = rep_vector(0, owt); // expected observations at each site in log scale
  vector[owt] exp_obs_log_v = rep_vector(0, owt); // expected observations at each site with modifier in log scale
  vector[n_ww_lab_sites] ww_site_mod; // site specific WW mod
  row_vector [ot + uot + ht] model_net_i; // number of net infected individuals shedding on each day (sum of individuals in dift stages of infection)
  real<lower=0> phi_h = inv_square(inv_sqrt_phi_h);
  vector<lower=0>[n_ww_lab_sites] sigma_ww_site;
  vector[n_weeks] log_r_mu_t_in_weeks; // log of state level mean R(t) in weeks
  vector<lower=0>[ot + ht] unadj_r; // state level R(t) before damping
  matrix[n_subpops, ot+ht] r_site_t; // site_level R(t)
  row_vector[ot + ht] unadj_r_site_t; // site_level R(t) before damping
  row_vector[ot + uot + ht] new_i_site; // site level incident infections per capita
  real<lower=0> pop_fraction; // proportion of state population that the subpopulation represents
  vector[ot + uot + ht] state_inf_per_capita = rep_vector(0, uot + ot + ht); // state level incident infections per capita
  matrix[n_subpops, ot + ht] model_log_v_ot; // expected observed viral genomes/mL at all observed and forecasted times
  real<lower=0> g = pow(log10_g, 10); // Estimated genomes shed per infected individual
  vector<lower=0, upper=1>[n_subpops] i_first_obs_over_n_site;
  // per capita infection incidence at the first observed time
  vector[n_subpops] initial_exp_growth_rate_site;
     // site level unobserved period growth rate

  // Site spatial trans params--------------------------------------------------
  real phi = exp(log_phi);
  real sigma_generalized = exp(log_sigma_generalized);
  real scaling_factor = exp(log_scaling_factor);
  matrix[n_subpops-1,n_subpops-1] non_norm_omega;
  matrix[n_subpops-1,n_subpops-1] norm_omega;
  matrix[n_subpops-1,n_subpops-1] sigma_mat;
  matrix[n_subpops-1,n_weeks] spatial_dev_ns_mat;
  matrix[n_subpops-1,n_weeks] log_r_site_t_in_weeks_matrix;
  vector[n_weeks] log_r_aux_site_t_in_weeks;
  matrix[n_subpops, n_weeks] combined_log_r_site_t_in_weeks;
  vector[n_weeks] log_r_site_t_in_weeks_vector;
  //----------------------------------------------------------------------------

  // State-leve R(t) AR + RW implementation:
  log_r_mu_t_in_weeks = diff_ar1(log_r_mu_intercept,
                                 autoreg_rt,
				 eta_sd,
				 w,
				 0);
  unadj_r = ind_m*log_r_mu_t_in_weeks;
  unadj_r = exp(unadj_r);

  // Shedding kinetics trajectory
  s = get_vl_trajectory(t_peak, viral_peak, dur_shed, gt_max);


  // Site level spatial Rt------------------------------------------------------
  if (corr_structure_switch == 0){
    // If no dist matrix given, use n_sites + 1 = n_subpops were all ind.
    non_norm_omega = independence_corr_func(n_subpops - 1);
    norm_omega = non_norm_omega;
  }
  else if (corr_structure_switch == 1){
    non_norm_omega = exponential_decay_corr_func(norm_dist_matrix, phi, l);
    norm_omega = matrix_normalization(non_norm_omega);
  }
  else if (corr_structure_switch == 2) {
    non_norm_omega = multiply_lower_tri_self_transpose(L_Omega);
    norm_omega = matrix_normalization(non_norm_omega);
  }
  else {
    reject("Model should not reach this point. Invalid corr_structure_switch value. Check model code");
  }
  sigma_mat = pow(sigma_generalized, 1.0 / (n_subpops - 1)) * norm_omega;
  for (i in 1:n_weeks) {
    spatial_dev_ns_mat[,i] = cholesky_decompose(sigma_mat) * non_cent_spatial_dev_ns_mat[,i];
  }
  log_r_site_t_in_weeks_matrix = construct_spatial_rt(
    log_r_mu_t_in_weeks,
    autoreg_rt_site,
    spatial_dev_ns_mat
  );
  //----------------------------------------------------------------------------
  // AUX site Rt----------------------------------------------------------------
  log_r_aux_site_t_in_weeks = construct_aux_rt(
    log_r_mu_t_in_weeks,
    autoreg_rt_site,
    scaling_factor,
    sqrt(sigma_generalized),
    norm_vec_aux_site,
    0
  );
  //----------------------------------------------------------------------------
  // Site Comb with AUX---------------------------------------------------------
  combined_log_r_site_t_in_weeks = append_row(log_r_site_t_in_weeks_matrix, log_r_aux_site_t_in_weeks');
  //----------------------------------------------------------------------------
  // Site level disease dynamics
  i_first_obs_over_n_site = inv_logit(logit(i_first_obs_over_n) +
      sigma_i_first_obs * eta_i_first_obs);
  initial_exp_growth_rate_site = mean_initial_exp_growth_rate +
     sigma_initial_exp_growth_rate * eta_initial_exp_growth_rate;

  for (i in 1:n_subpops) {
    real log_i0_site = log(i_first_obs_over_n_site[i]) - uot * initial_exp_growth_rate_site[i];
     //convert from weekly to daily
     log_r_site_t_in_weeks_vector = to_vector(combined_log_r_site_t_in_weeks[i, :]);
     unadj_r_site_t = exp(to_row_vector(ind_m*(log_r_site_t_in_weeks_vector)));

    {
      tuple(vector[num_elements(state_inf_per_capita)], vector[num_elements(unadj_r)]) output;
      output = generate_infections(
        to_vector(unadj_r_site_t),
	      uot,
	      gt_rev_pmf,
	      log_i0_site ,
	      initial_exp_growth_rate_site[i],
	      ht,
        infection_feedback,
	      infection_feedback_rev_pmf
      );
      new_i_site = to_row_vector(output.1);
      r_site_t[i] =  to_row_vector(output.2);
    }

    // For each site, tack on number of state infections
    // site level infection dynamics sum to the total state infections:
    pop_fraction = subpop_size[i] / norm_pop;
    state_inf_per_capita +=  pop_fraction * to_vector(new_i_site);

    model_net_i = to_row_vector(convolve_dot_product(to_vector(new_i_site),
                               reverse(s), (uot + ot + ht)));


    model_log_v_ot[i] = log(10) * log10_g +
      log(model_net_i[(uot+1):(uot + ot + ht) ] + 1e-8) -
      log(mwpd);
  }


  // Set up p_hosp as an AR(1) process that regresses back towards the initial value of p_hosp
	p_hosp = assemble_p_hosp(p_hosp_m, p_hosp_mean, p_hosp_w_sd,
                           autoreg_p_hosp, p_hosp_w, tot_weeks, 1);

  // Expected hospital admissions per capita:
  // This is a convolution of incident infections and the hospital-admission delay distribution
  // generates all hospitalizations, across unobserved time, observed time, and forecast time
  model_hosp_per_capita = convolve_dot_product(p_hosp .* state_inf_per_capita, reverse(inf_to_hosp),
                                    ot + uot + ht);

  // predicted hospital admissions per capita
  exp_obs_hosp_per_capita_no_wday_effect = model_hosp_per_capita[uot + 1 : uot + ot];
  // apply the weekday effect so these are distributed with fewer admits on Sat & Sun
  // multiply by state population to convert from predicted per capita admissions to
  // predicted absolute admissions
  exp_obs_hosp = state_pop * day_of_week_effect(
  	       exp_obs_hosp_per_capita_no_wday_effect[hosp_times],
	       day_of_week[hosp_times],
	       hosp_wday_effect);

  // Observations at the site level (genomes/person/day) are:
  // get a vector of genomes/person/day on the days WW was measured
  // These are the true expected genomes at the site level before observation error
  // (which is at the lab-site level)
  for (i in 1:owt) {
    exp_obs_log_v_true[i] = model_log_v_ot[ww_sampled_sites[i], ww_sampled_times[i]];
  }

  // modify by lab-site specific variation (multiplier!)
  ww_site_mod = ww_site_mod_raw * ww_site_mod_sd;
  // LHS log transformed obs genomes per person-day, RHS multiplies the expected observed
  // genomes by the site-specific multiplier at that sampling time
  exp_obs_log_v = exp_obs_log_v_true + ww_site_mod[ww_sampled_lab_sites];
  // Option to add a population offset here at some point  log(model_V) + site_level_multiplier+ pop_ww[ww_sampled_sites]

  // Get the transformed lab-site level error
  // log(sigma_site) ~ Normal(log(mode_sigma_site), sd_log_sigma_site)
  sigma_ww_site = exp(
  log(mode_sigma_ww_site) + sd_log_sigma_ww_site * eta_log_sigma_ww_site);

}

// Prior and sampling distribution
model {
  // priors
    // for spatial--------------------------------------------------------------
    to_vector(non_cent_spatial_dev_ns_mat) ~ std_normal();
    norm_vec_aux_site ~ std_normal();
    log_sigma_generalized ~ normal(log_sigma_generalized_mu_prior, log_sigma_generalized_sd_prior);
    log_phi ~ normal(log_phi_mu_prior, log_phi_sd_prior);
    log_scaling_factor ~ normal(log_scaling_factor_mu_prior, log_scaling_factor_sd_prior);
    if (corr_structure_switch == 2){
      L_Omega ~ lkj_corr_cholesky(0.25);
    }
    //--------------------------------------------------------------------------

  w ~ std_normal();
  eta_sd ~ normal(0, eta_sd_sd);
  autoreg_rt_site ~ beta(autoreg_rt_site_a, autoreg_rt_site_b);

  autoreg_rt ~ beta(autoreg_rt_a, autoreg_rt_b);
  autoreg_p_hosp ~ beta(autoreg_p_hosp_a, autoreg_p_hosp_b);
  log_r_mu_intercept ~ normal(r_logmean, r_logsd);
  sigma_rt ~ normal(0, sigma_rt_prior);
  i_first_obs_over_n ~ beta(i_first_obs_over_n_prior_a,
                             i_first_obs_over_n_prior_b);
  sigma_i_first_obs ~ normal(sigma_i_first_obs_prior_mode,
		             sigma_i_first_obs_prior_sd);
  eta_i_first_obs ~ std_normal();
  sigma_initial_exp_growth_rate ~ normal(sigma_initial_exp_growth_rate_prior_mode,
                                         sigma_initial_exp_growth_rate_prior_sd);
  eta_initial_exp_growth_rate ~ std_normal();
  mean_initial_exp_growth_rate ~ normal(mean_initial_exp_growth_rate_prior_mean, mean_initial_exp_growth_rate_prior_sd);
  inv_sqrt_phi_h ~ normal(inv_sqrt_phi_prior_mean, inv_sqrt_phi_prior_sd);
  mode_sigma_ww_site ~ normal(mode_sigma_ww_site_prior_mode,
                              mode_sigma_ww_site_prior_sd);
  sd_log_sigma_ww_site ~ normal(sd_log_sigma_ww_site_prior_mode,
                                sd_log_sigma_ww_site_prior_sd);
  eta_log_sigma_ww_site ~ std_normal();
  log10_g ~ normal(log10_g_prior_mean, log10_g_prior_sd);
  hosp_wday_effect ~ dirichlet(hosp_wday_effect_prior_alpha);
  p_hosp_mean ~ normal(logit(p_hosp_prior_mean), p_hosp_sd_logit); // logit scale
  p_hosp_w ~ std_normal();
  p_hosp_w_sd ~ normal(0, p_hosp_w_sd_sd);
  t_peak ~ normal(t_peak_mean, t_peak_sd);
  viral_peak ~ normal(viral_peak_mean, viral_peak_sd);
  dur_shed ~ normal(dur_shed_mean, dur_shed_sd);
  ww_site_mod_raw ~ std_normal();
  ww_site_mod_sd ~ normal(0, ww_site_mod_sd_sd);
  infection_feedback ~ lognormal(inf_feedback_prior_logmean, inf_feedback_prior_logsd);

  //Compute log likelihood
  if (compute_likelihood == 1) {
    if (include_ww == 1) {
      // Both genomes/person/day and observation error are now vectors
      //log_conc ~ normal(exp_obs_log_v, sigma_ww_site[ww_sampled_lab_sites]);
      // if non-censored: P(log_conc | expected log_conc)
      log_conc[ww_uncensored] ~ normal(exp_obs_log_v[ww_uncensored], sigma_ww_site[ww_sampled_lab_sites[ww_uncensored]]);
      // The stdev is at the lab-site-level
      // if censored: P(expected_log_conc <= LOD)
      target +=  normal_lcdf(ww_log_lod[ww_censored]| exp_obs_log_v[ww_censored],
                            sigma_ww_site[ww_sampled_lab_sites[ww_censored]]);
    }

    if (include_hosp == 1) {
      hosp ~ neg_binomial_2(exp_obs_hosp, phi_h);
    }
  } // end if for computing log likelihood
}

generated quantities {
  array[ot + ht] real pred_hosp;
  array[ot + ht] real pred_new_i;
  array[n_ww_lab_sites, ot + ht] real pred_ww; /// viral genome copies/person/day
  vector[ot + ht] exp_state_ww_conc;
  vector[ot + ht] state_log_c;
  vector[uot + ot + ht] state_model_net_i;
  vector<lower=0>[ot + ht] rt; // state level R(t)

  pred_hosp = neg_binomial_2_rng(state_pop * day_of_week_effect(model_hosp_per_capita[uot + 1 :
                                                        uot + ot + ht],
                                                        day_of_week,
                                                        hosp_wday_effect),
                                 phi_h);
  pred_new_i = neg_binomial_2_rng(state_pop * state_inf_per_capita[uot + 1 : uot + ot + ht], phi_h);

  // Here need to iterate through each lab-site, find the corresponding site
  // and apply the expected lab-site error
  for(i in 1:n_ww_lab_sites) {
    pred_ww[i] = normal_rng(model_log_v_ot[lab_site_to_site_map[i], 1 : ot + ht] + ww_site_mod[i],
                            sigma_ww_site[i]);
  }

  state_model_net_i = convolve_dot_product(state_inf_per_capita,
                                           reverse(s), (uot + ot + ht));
  state_log_c = log(10) * log10_g +
    log(state_model_net_i[(uot + 1): (uot + ot + ht)] + 1e-8) -
    log(mwpd);

  exp_state_ww_conc = exp(state_log_c);

  // Deterministic calculation of state level R(t) from incident infections
  // and the generation interval, with a 0 added to account for the fact
  // that our forward process assumes no contribution to onwards infections
  // on the day of infection
  rt = (state_inf_per_capita ./ convolve_dot_product(state_inf_per_capita,
    reverse(append_row(0,generation_interval)), uot + ot + ht))[uot+1: uot + ot + ht];

}
