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
  real<lower=0> a;
  real<lower=0> b;
  int<lower = 0> D;
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
  target += uniform_lpdf(u | -5, 5);
  target += -((u[D]^2 + u[1]^2 - a)^2 / b)^2 - 
    sum(square((square(u[1:(D-1)]) + square(u[2:D]) - a) / b));
}

