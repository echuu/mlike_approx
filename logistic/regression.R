


library(rstan)
library(rstanarm)
library(ggplot2)
library(bayesplot)
library(dplyr)

options(mc.cores = parallel::detectCores()) 

# http://biostat.mc.vanderbilt.edu/wiki/pub/Main/StatisticalComputingSeries/bayes_reg_rstanarm.html


data(cars)

# specify probability model (in this case a linear regression model)
# rstanarm has default priors in case they are not user-specified
# stan_glm() used to draw samples from the posterior
# stan_trace() fo trace plots, which show sequential draws from the posterior

# model     : dist = beta0 + beta1 * speed + eps, eps ~ N(0, sigma^2)
# parameters: beta0, beta1, sigma^2
glm_post1 = stan_glm(dist ~ speed, data = cars, family = gaussian, refresh = 0)

# look at trace plots for each of the parameters
stan_trace(glm_post1, pars=c("(Intercept)","speed","sigma"))

summary(glm_post1)

# can also look at PP: if model is a good fit, then we should be able to use
# it to generate data that looks a lot like the data we observed
pp_check(glm_post1)


# another way to look at posterior predictive checks
ppc_intervals(
 y = cars$dist,
 yrep = posterior_predict(glm_post1),
 x = cars$speed)

# change the choice in prior using `prior` argument
glm_post2 = stan_glm(dist ~ speed, data = cars, family = gaussian,
                     prior = normal(2, 0.5, autoscale = FALSE))

summary(glm_post2)

# obtain samples from the posterior
as.data.frame(glm_post1) %>% head
as.data.frame(glm_post2) %>% head


# ------------------------------------------------------------------------------

n = 500; p = 50;
X = matrix(rnorm(n * p), nrow = 500)
Eps = rnorm(500, 0, 5)
beta = rep(4.5, p)
Y = X %*% beta + Eps;
dat = data.frame(Y = Y, X = X)
fit_lin <- stan_glm(Y ~ X, data = dat, family = gaussian())
summary(fit_li)n
post_samps = as.data.frame(fit_lin)



# bayesian logistic regression model -------------------------------------------

set.seed(1)
N = 300
D = 3
X = matrix(rnorm(N * D), nrow = N)
beta = sample(c(-5:5), D, replace = TRUE) # (-1, 5, 2) for D = 3
z = X %*% beta
pr = exp(z) / (1 + exp(z))
y = rbinom(N, 1, pr)

df = data.frame(y, X)

glm_fit = glm(y ~ . -1, data = df, family = 'binomial')
summary(glm_fit)


# fit bayesian version with N(0, tau^(-1) I_D) prior on beta

bglm_fit = stan_glm(y ~ . -1, data = df, 
                    prior = normal(location = 0, scale = 4),
                    family = binomial(link = 'logit'))
summary(bglm_fit)

as.data.frame(bglm_fit) %>% head


stan_trace(bglm_fit, pars=c("X1","X2","X3"))









