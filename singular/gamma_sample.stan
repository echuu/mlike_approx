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


data {
  int<lower=0> N; // sample size
}


parameters {
  vector<lower=0, upper=1>[2] u;   // 2-dim parameter
}


model {
  
  target += uniform_lpdf(u[1] | 0, 1) + uniform_lpdf(u[2] | 0, 1); // log prior
  target += -N * square(u[1]) * square(square(u[2])); // log likelihood
  
}

