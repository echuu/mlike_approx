

// The input data is a vector 'y' of length 'N'.
data {
  // int<lower=0> N;
  int<lower=0> p;
  
  real a_n;
  real b_n;
  vector[p] mu_star;
  matrix[p, p] V_star;
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[p] beta;
  real<lower=0> sigmasq;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // y ~ normal(mu, sigma);
  target += inv_gamma_lpdf(sigmasq | a_n, b_n);
  target += multi_normal_lpdf(beta | mu_star, sigmasq * V_star);
  
}

