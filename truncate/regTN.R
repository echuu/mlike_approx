

library(TruncatedNormal)
library(tmg)
library(mvtnorm)

LEN_PATH  = "C:/Users/ericc/mlike_approx"
setwd(LEN_PATH)

source("partition/partition.R")
source("extractPartition.R")
source("hybrid_approx.R")
source("truncate/regTN_helper.R")
D = 2
N = 200
I_D = diag(1, D)


# prior mean
mu_0 = rep(0, D)

tau     = 1 / 4          # precision: inverse of variance
sigmasq = 4              # true variance (1 x 1) 

# true value of beta
set.seed(1)
beta = sample(0:10, D, replace = T)
beta = runif(D)
beta = c(runif(D-1, 0, 1), 0)

# generate the regression data -------------------------------------------------

X   = matrix(rnorm(N * D), N, D)                # (N x D) design matrix
eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))    # (N x 1) errors vector
y   = X %*% beta + eps                          # (N x 1) response vector

# compute posterior parameters -------------------------------------------------
Q_beta     =  1 / sigmasq * (t(X) %*% X + tau * diag(1, D))
Q_beta_inv =  solve(Q_beta)
b          =  1 / sigmasq * t(X) %*% y
mu_beta    =  Q_beta_inv %*% b


# create prior, post objects to be passed into the hml algorithm 
prior = list(y = y, X = X, sigmasq = sigmasq, tau = tau, N = N, D = D,
             Q_beta = Q_beta, b = b)
post = list(Q_beta = Q_beta, Q_beta_inv = Q_beta_inv, mu_beta = mu_beta, b = b)


# sample from posterior distribution -------------------------------------------

## TODO: udpate the parameters that are passed into rtmg() function
## i think we need the **precision** matrix, not the covariance matrix

# R       =  rep(0, D)   # linear part R'X
# ff      =  diag(D)     # ?
# gg      =  rep(0, D)   # constraints (truncate to positive orthant)
# initial =  rep(1, D)   # Set initial point for the Markov chain


J        = 100         # number of MC samples
B        = 10          # number of batch estimates
N_approx = 1           # number of estimates to compute per iteration b

# samples = data.frame(rtmg(J * B,           # number of samples
#                           M = Q_beta_inv,  # (D x D) precision matrix
#                           r = R,           # linear coeff of log-density
#                           initial,         # initial value of markov chain
#                           f = ff,          #
#                           g = gg))         # m-dim with constraints

# plot 2-d samples to verify truncation
# plot(samples)

samples = data.frame(rtmvnorm(J * B, c(mu_beta), Q_beta_inv, rep(0, D), rep(Inf, D)))
# samples = data.frame(rtmvnorm(J * B, R, Q_beta_inv, rep(0, D), rep(Inf, D)))

plot(samples)

## TODO: be able to evaluate log-prior of each of the posterior samples
# logprior = TruncatedNormal::dtmvnorm(samples[1,], mu = c(mu_beta),
#                                      sigma = Q_beta_inv, lb = rep(0, D),
#                                      ub = rep(Inf, D), log = T)

# dtmvnorm(unlist(unname(samples[5,])), mu = c(mu_beta),
#          sigma = Q_beta_inv, lb = rep(0, D),
#          ub = rep(Inf, D), log = T)
# 
# dmvnorm(unlist(unname(samples[5,])), mean = c(mu_beta),
#         sigma = Q_beta_inv, log = T)
# 
# 
# ## compute probability to compare to closed form
# 
# TruncatedNormal::pmvnorm(c(mu_beta), Q_beta_inv, lb = rep(-Inf, D), ub = rep(Inf, D))[1]  
# TruncatedNormal::pmvnorm(c(mu_beta), Q_beta_inv, lb = rep(0, D), ub = rep(Inf, D))[1]  
# 
# u = unlist(unname(samples[5,]))
# dtmvnorm(u, mu = rep(0, D),
#          sigma = tau / sigmasq * I_D , lb = rep(0, D),
#          ub = rep(Inf, D), log = T)
# 
# dtmvnorm(u, mu = rep(0, D),
#          sigma = tau / sigmasq * I_D , lb = rep(0, D),
#          ub = rep(Inf, D))

(lil_0 = -0.5 * N * log(2 * pi) - 0.5 * (N + D) * log(sigmasq) + 
    0.5 * D * log(tau) - 0.5 * log_det(Q_beta) - 
    1 / (2 * sigmasq) * sum(y^2) + 0.5 * sum(b * mu_beta))


# lil_0 - log(TruncatedNormal::pmvnorm(rep(0, D), tau / sigmasq * I_D, lb = rep(0, D), ub = rep(Inf, D))[1]) 

lil_0 + log(TruncatedNormal::pmvnorm(mu_beta, Q_beta_inv, lb = rep(0, D), ub = rep(Inf, D))[1]) 


# psi(u, prior)
# 
# - sum(dnorm(y, X %*% u, sqrt(sigmasq), log = T)) - 
#     dtmvnorm(u, mu = rep(0, D), sigma = tau / sigmasq * diag(1, D) , 
#              lb = rep(0, D), ub = rep(Inf, D), log = T)
# 

# compute *true* log marginal likelihood ---------------------------------------
# prob = TruncatedNormal::pmvnorm(mu, Sigma, lb = rep(0, D), ub = rep(Inf, D))     
# prob[1] %>% log

# run algorithm ----------------------------------------------------------------


set.seed(1)


TruncatedNormal::dtmvnorm(c(0.8210398, 0.1810452), mu = rep(0, D),
                 sigma = sigmasq / tau * diag(1, D) , lb = rep(0, D),
                 ub = rep(Inf, D), log = F)

prod(dnorm(c(0.8210398, 0.1810452), 0, sqrt(sigmasq / tau), log = F)) * 4


# start_time <- Sys.time()
u_df = preprocess(samples, D, prior)

psi_fast = u_df$psi_u
head(psi_fast)


sum(psi_fast != psi_slow)
psi_fast[psi_fast != psi_slow]
which(psi_fast != psi_slow)

samples[35,]

psi_slow = u_df$psi_u
head(psi_slow)




# end_time <- Sys.time()
# end_time - start_time

# u_df %>% head

hml_approx = hml(N_approx, D, u_df, 100, prior)
hml_approx$hybrid_vec













