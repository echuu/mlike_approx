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
  real<lower=0> mu;
  real<lower=0> sigma;
  int<lower=0> D;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[D] u;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  target += uniform_lpdf(u | -100, 100);
   //target += D * log(0.5) + sum(log(exp(cauchy_lpdf(u | mu, sigma)) + exp(cauchy_lpdf(u | -mu, sigma))));
  for (d in 1:D)
    target += log(0.5) + log(exp(cauchy_lpdf(u[d] | mu, sigma)) + exp(cauchy_lpdf(u[d] | -mu, sigma)));
}

