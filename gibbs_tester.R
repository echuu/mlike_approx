

p = 8             # number of columns in X
q = 6             # number of columns in Y
r = 2             # number of columns in B and A
D = r * p + q * r # dimension of each MCMC sample
n = 100           # number of rows in X and Y
sig2 = 10^(-2)          # fixed for now.
del = 10^(-2)     # prior parameter -- one of these is squared version ?

set.seed(1)
A_0 = matrix(rnorm(p * r, 0, 1), p, r) # (p x r) matrix
B_0 = matrix(rnorm(q * r, 0, 1), q, r) # (q x r) matrix

nMCMC = 300       # number of MCMC samples from the posterior AFTER burnin
nBurn = 500       # number of samples to discard

# sample from the posterior
gibbs_obj = sampleRRR(nMCMC, nBurn, A_0, B_0, p, q, r, r, D, n, sig2, del)

# extract posterior samples from gibbs object
u_samps = gibbs_obj$u_samps


# u = u_samps[300,] %>%  unname %>% unlist
# 
# matrix(u[1:(p * r)], p, r)
# t(matrix(u[(p * r + 1):D], r, q))
# 
# A_0 %*% t(B_0)
# 
# matrix(u[1:(p * r)], p, r) %*% matrix(u[(p * r + 1):D], r, q)

param_list = list(p = p, q = q, r = r, n = n, d = D,  # dimensions variables
                  Y = gibbs_obj$Y, X = gibbs_obj$X,   # response, design matrix
                  XtX = gibbs_obj$XtX, Xty = gibbs_obj$Xty,
                  sig2 = sig2, del = del)             # prior params

# extract posterior samples from gibbs object
# u_samps = gibbs_obj$u_samps
# 
# 
# u = u_samps[300,] %>%  unname %>% unlist
# 
# matrix(u[1:(p * r)], p, r)
# matrix(u[(p * r + 1):D], r, q)
# 
# 
# u0 = u_samps[1,] %>% unname %>% unlist()
# rrr_logprior(u0, param_list)
# rrr_loglik(u0, param_list)


# evaluate psi(u) for each of the posterior samples
u_df = preprocess(u_samps, D, param_list) # J x (d + 1) 


ll_max = loglik_true(A_0, B_0, param_list)

# generate hybrid approximation
hml_approx = hml(1, D, u_df, nMCMC, param_list)
hml_approx$hybrid_vec - ll_max

hml_approx$n_taylor





