


J = 500                  # num of MCMC samples (from true posterior in this case)
N = 500                  # number of observations
D = 5                    # num rows/cols in the covariance matrix
D_u = 0.5 * D * (D + 1)  # dimension of u that is fed into the tree


## wishart prior parameters
Omega = diag(1, D)  # scale matrix
nu    = D + 1       # degrees of freedom


set.seed(1)
Sigma = matrix(rWishart(1, D, Omega), D)
is.positive.definite(Sigma)


# obtain lower cholesky factor
L = t(chol(Sigma))
u = L[lower.tri(L, diag = T)]

# reconstruct L
Lp = matrix(0, D, D)
Lp[lower.tri(L, diag = T)] = u


# (1) use cov_logprior_sigma() to evalute Sigma
cov_logprior_sigma(Sigma, param_list)

# (2) log(diwishart(LL')) to evaluate Sigma
log(diwish(L %*% t(L), nu, Omega))


# (1) use cov_logprior() to evaluate u
cov_logprior(u, param_list)


# (2) use log(diwishart(LL')) + Jacobian term
logDiagL = log(diag(Lp))         # log of diagonal terms of L
logJacTerm = D * log(2) + sum((D + 1 - 1:D) * logDiagL)

log(diwish(Lp %*% t(Lp), nu, Omega)) + logJacTerm







# (1) and (2) should match



















