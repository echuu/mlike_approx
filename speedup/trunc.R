
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")

# look in this file for psi()
source("C:/Users/ericc/mlike_approx/truncate/regTN_helper.R")

set.seed(123)
D = 5
N = 200
I_D = diag(1, D)

n_samps = 10
J       = 1000
B       = 1000 # number of replications


#### generate data -------------------------------------------------------------

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


samples = rtmvnorm(J, c(mu_beta), Q_beta_inv, rep(0, D), rep(Inf, D))
u_df = preprocess(data.frame(samples), D, prior)


head(u_df)

hml_approx = hml_const(1, D, u_df, J, prior)
hml_approx$const_vec ## check this output after doing slow vs. fast psi




