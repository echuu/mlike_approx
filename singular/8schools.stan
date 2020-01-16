//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> J;            // number of schools
  real y[J];                 // estimated treatment effects
  real<lower=0> sigma[J];    // standard error of effect estimates 
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {    
  real mu;                   // population treatment effects
  real<lower=0> tau;         // standard deviation in treatment effects
  vector[J] eta;             // unscaled deviation from mu by school
}

transformed parameters {
  vector[J] theta = mu + tau * eta;    // school treatment effects
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // y ~ normal(mu, sigma);
  
  target += normal_lpdf(eta | 0, 1);        // prior log-density
  target += normal_lpdf(y | theta, sigma);  // log-likelihood
  
}

