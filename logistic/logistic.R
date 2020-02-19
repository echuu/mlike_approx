
library(rstan)
library(rstanarm)
library(ggplot2)
library(bayesplot)
library(dplyr)
library(gtools)    # for inv.logit() function

options(mc.cores = parallel::detectCores()) 

# save compiled stan program to hard disk so that don't need to recompile
rstan_options(auto_write = TRUE)

LEN_PATH  = "C:/Users/ericc/mlike_approx"
setwd(LEN_PATH)


source("partition/partition.R")
source("extractPartition.R")
source("hybrid_approx.R")
source("logistic/logistic_helper.R")


set.seed(1)

N = 300
D = 3
X = matrix(rnorm(N * D), nrow = N)
beta = sample(c(-5:5), D, replace = TRUE) # (0, 1, 3) for D = 3
z = X %*% beta
pr = exp(z) / (1 + exp(z))
y = rbinom(N, 1, pr)

df = data.frame(y, X)

# prior precision for beta
tau = 0.5
prior = list(y = y, X = X, D = D, N = N, tau = tau)

# number of MCMC samples
J = 4000

## TODO: figure out how to specify the number of posterior samples used
bglm_fit = stan_glm(y ~ . -1, data = df, 
                    prior = normal(location = 0, scale = 1/sqrt(tau)),
                    family = binomial(link = 'logit'),
                    refresh = 0)
summary(bglm_fit)

u_samps = bglm_fit %>% as.data.frame() # extract posterior samples

psi_u = apply(u_samps, 1, psi, prior = prior)
u_df = preprocess(u_samps, D, prior)

hml_approx = hml(1, D, u_df, J, prior) 

hml_approx$const_vec
hml_approx$taylor_vec
hml_approx$hybrid_vec

hml_approx$partition


# LIL_const[,k]  = hml_approx$const_vec
# LIL_taylor[,k] = hml_approx$taylor_vec
# LIL_hybrid[,k] = hml_approx$hybrid_vec





