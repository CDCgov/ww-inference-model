library(tidyverse)
library(rstan)



stan_model_code_non_cent <- "
functions {
  matrix exponential_decay_corr_func(matrix dist_matrix, real phi, real l) {
    int k = rows(dist_matrix);
    matrix[k,k] corr_matrx = exp(-pow(dist_matrix / phi, l));
    return corr_matrx;
  }
}

data {
  int<lower=1> n_sites;
  int<lower=1> n_times;
  matrix[n_sites,n_sites] dist_matrix;
  real<lower=0> l;
  matrix[n_sites, n_times] obs_data;
}

parameters {
  matrix[n_sites,n_times] alpha;
  real<lower=0> phi;
  real<lower=0> sigma_eps;
  //real<lower=0> sigma_obs; Can get rid of this and assume iid normal
}

transformed parameters {
  matrix[n_sites,n_sites] Sigma;
  matrix[n_sites,n_times] epsilon;

  Sigma = sigma_eps^2 * exponential_decay_corr_func(dist_matrix, phi, l);

  for (i in 1:n_times) {
    epsilon[,i] = cholesky_decompose(Sigma) * alpha[, i];
  }
}

model {
  //Priors
  to_vector(alpha) ~ std_normal();
  phi ~ uniform(0, 50);
  sigma_eps ~ gamma(sqrt(0.1),1);

  // Likelihood
  // Assumes normally distributed observation error.
  // obs_data are your vectors of draws the MVN at each time point (column
  // bound to produce a matrix). The epsilon matrix are the
  // mean of the expected realizations of the MVN draws.
  for (i in 1:n_times) {
    obs_data[,i] ~ normal(epsilon[,i], 1);
  }
}
"


stan_model_code_cent <- "
functions {
  matrix exponential_decay_corr_func(matrix dist_matrix, real phi, real l) {
    int k = rows(dist_matrix);
    matrix[k,k] corr_matrx = exp(-pow(dist_matrix / phi, l));
    return corr_matrx;
  }
}

data {
  int<lower=1> n_sites;
  int<lower=1> n_times;
  matrix[n_sites,n_sites] dist_matrix;
  real<lower=0> l;
  matrix[n_sites,n_times] obs_data;
}

parameters {
  real<lower=0> phi;
  real<lower=0> sigma_eps;
}

transformed parameters {
  matrix[n_sites,n_sites] Sigma;

  Sigma = sigma_eps^2 * exponential_decay_corr_func(dist_matrix, phi, l);
}

model {
  //Priors
  phi ~ uniform(0, 50);
  sigma_eps ~ gamma(sqrt(0.1),1);

  // Likelihood
  for (i in 1:n_times) {
    obs_data[,i] ~ multi_normal(rep_vector(0, n_sites), Sigma);
  }
}
"

stan_model_non_cent <- stan_model(
  model_code = stan_model_code_non_cent
)
stan_model_cent <- stan_model(
  model_code = stan_model_code_cent
)

# data
n_sites <- 4
n_times <- 12
dist_matrix <- as.matrix(
  dist(
    data.frame(
      x = c(85, 37, 48, 7),
      y = c(12, 75, 81, 96),
      diag = TRUE,
      upper = TRUE
    )
  )
)
l <- 1

# obs_data generation
phi_true <- 25
sigma_eps_true <- sqrt(0.1)
sigma_matrix_true <- sigma_eps_true^2 * exponential_decay_corr_func_r(list(
  dist_matrix = dist_matrix,
  phi = phi_true,
  l = l
))
# Generate samples from a MVM with mean 0 and covariance matrix sigma
# epsilon_t ~ MVN(0, sigma_matrix) # nolint
# for each time point. These are your observations!
obs_data <- matrix(data = 0, nrow = n_sites, ncol = n_times)
for (i in 1:n_times) {
  obs_data[, i] <- t(mvrnorm(
    n = 1,
    mu = rep(0, n_sites),
    Sigma = sigma_matrix_true
  ))
}

# Save the data in a list for Stan
data_list <- list(
  n_sites = n_sites,
  n_times = n_times,
  dist_matrix = dist_matrix,
  l = l,
  obs_data = obs_data
)

# Fit the model
fit_non_cent <- sampling(
  stan_model_non_cent,
  data = data_list,
  iter = 5000,
  chains = 1
) # 2.569sec
fit_cent <- sampling(
  stan_model_cent,
  data = data_list,
  iter = 5000,
  chains = 1
) # 3.083sec


# averaging eps time points
fit_non_cent_extract <- rstan::extract(fit_non_cent)$epsilon
epsilon_non_cent <- data.frame()
fit_cent_extract <- rstan::extract(fit_cent)$epsilon
epsilon_cent <- data.frame()
for (i in 1:n_times) {
  temp <- fit_non_cent_extract[, , i] %>%
    as.data.frame() %>%
    summarise(across(everything(), ~ mean(.x)))
  epsilon_non_cent <- rbind(epsilon_non_cent, temp)
  temp <- fit_cent_extract[, , i] %>%
    as.data.frame() %>%
    summarise(across(everything(), ~ mean(.x)))
  epsilon_cent <- rbind(epsilon_cent, temp)
}


epsilon_non_cent <- epsilon_non_cent %>%
  as.matrix() %>%
  t()
epsilon_cent <- epsilon_cent %>%
  as.matrix() %>%
  t()

cramer::cramer.test(
  epsilon_non_cent,
  epsilon_cent
)
