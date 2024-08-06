library(tidyverse)
library(MASS)

# For this example, we will be exploring noncentering vs centering a MVN
# using the Choleskey decomposition functionality in Stan.
# The centered model we will be using is the following :
# // epsilon ~ MVN(0, Sigma)
# // obs_data ~ MVN(epsilon, I)
# The non-centered model will be using the following :
# // alpha ~ N(0,1)
# // Sigma = LL'
# // epsilon = 0 + L * alpha
# // obs_data ~ MVN(epsilon, I)


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
  matrix[n_sites, n_times] epsilon;
}

transformed parameters {
  matrix[n_sites,n_sites] Sigma;

  Sigma = sigma_eps^2 * exponential_decay_corr_func(dist_matrix, phi, l);
}

model {
  //Priors
  phi ~ uniform(0, 50);
  sigma_eps ~ gamma(sqrt(0.1),1);
  for (i in 1:n_times) {
    epsilon[,i] ~ multi_normal(rep_vector(0, n_sites), Sigma);
  }

  // Likelihood
  for (i in 1:n_times) {
    obs_data[,i] ~ normal(epsilon[,i], 1);
  }
}
"

stan_model_non_cent <- cmdstanr::write_stan_file(
  stan_model_code_non_cent
)
stan_model_cent <- cmdstanr::write_stan_file(
  stan_model_code_cent
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
epsilon <- matrix(data = 0, nrow = n_sites, ncol = n_times)
obs_data <- matrix(data = 0, nrow = n_sites, ncol = n_times)
for (i in 1:n_times) {
  epsilon[, i] <- t(mvrnorm(
    n = 1,
    mu = rep(0, n_sites),
    Sigma = sigma_matrix_true
  ))
  obs_data[, i] <- t(mvrnorm(
    n = 1,
    mu = epsilon[, i],
    Sigma = diag(1, n_sites, n_sites)
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
model_non_cent <- cmdstanr::cmdstan_model(
  stan_model_non_cent
)
fit_non_cent <- model_non_cent$sample(
  data_list
)
model_cent <- cmdstanr::cmdstan_model(
  stan_model_cent
)
fit_cent <- model_cent$sample(
  data_list
)


# averaging eps time points

fit_non_cent_extract <- fit_non_cent$draws(variables = "epsilon")
fit_cent_extract <- fit_cent$draws(variables = "epsilon")
epsilon_non_cent <- data.frame()
epsilon_cent <- data.frame()
for (i in 1:n_times) {
  temp_non_cent <- fit_non_cent_extract[, , i] %>%
    as.data.frame() %>%
    summarise(across(everything(), ~ mean(.x))) %>%
    `colnames<-`(c(
      "eps1", "eps2", "eps3", "eps4"
    ))
  epsilon_non_cent <- rbind(epsilon_non_cent, temp_non_cent)
  temp_cent <- fit_cent_extract[, , i] %>%
    as.data.frame() %>%
    summarise(across(everything(), ~ mean(.x))) %>%
    `colnames<-`(c(
      "eps1", "eps2", "eps3", "eps4"
    ))
  epsilon_cent <- rbind(epsilon_cent, temp_cent)
}


epsilon_non_cent <- epsilon_non_cent %>%
  as.matrix()
epsilon_cent <- epsilon_cent %>%
  as.matrix()

cramer::cramer.test(
  epsilon_non_cent,
  epsilon_cent
)
