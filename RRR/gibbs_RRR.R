


library(mvtnorm)
library(dplyr)




sampleRRR = function(n_samps, n_burn, A_0, B_0, p, q, r_0, r, D, N, sig2, del) {
    
    # identity matrices
    I_q = diag(1, q)
    I_r = diag(1, r)
    
    eps = matrix(rnorm(N * q, 0, sqrt(sig2)), N, q)   # (n x q) error
    X = rmvnorm(N, mean = rep(0, p), diag(1, p))      # (n x p) response
    Y   = X %*% A_0 %*% t(B_0) + eps                  # (n x q) response 
    # C   = A_0 %*% t(B_0)                            # (p x q)
    
    XtX = t(X) %*% X
    Xty = t(X) %*% Y
    
    B = matrix(rnorm(q * r, 0, 1), q, r) # (q x r) matrix starting point MCMC
    A = matrix(rnorm(p * r, 0, 1), p, r) # (p x r)
    
    # A = A_0 # (p x r) matrix for starting point of MCMC
    
    n_iters = n_samps + n_burn
    u_samps = matrix(0, n_iters, D)
    
    for (g in 1:n_iters) {
        
        ## (1)  sample from the conditional distribution: B' | - 
        Btvarpart = sig2 * solve(t(A) %*% XtX %*% A + del^2 * I_r)
        Btvar = I_q %x% Btvarpart
        Btmu = Btvarpart %*% t(A) %*% Xty / sig2
        
        Bt_row = c(rmvnorm(1, c(Btmu), Btvar)) # (1 x rq) row vector
        
        ## even though Bt is not stored in the posterior sample df in 
        ## matrix form, we still need Bt in matrix form to compute
        ## the conditional distribution of A | - in the next part
        Bt = matrix(Bt_row, r, q)              # (r x q) matrix
        B = t(Bt)
        
        Mt = solve(t(B) %*% B) %*% (t(B) %*% t(Y) %*% X) %*% solve(XtX)
        M = t(Mt)
        
        BtB = t(B) %*% B
        
        # ----------------------------------------------------------------------
        
        ## (2) sample from the conditional distribution A | -
        Avarpart = (BtB %x% XtX / sig2)
        Avar = solve(del^2 / sig2 * diag(1, nrow(Avarpart)) + Avarpart)
        Amu = Avar %*% Avarpart %*% c(M)
        
        A_row = c(rmvnorm(1, Amu, Avar))
        
        ## even though A is not stored in the posterior sample df in matrix form
        ## we still need A in matrix form to compute the conditional 
        ## distribution of B | - in the next iteration
        A = matrix(A_row, p, r) # (p x r) matrix
        
        u_samps[g,] = c(A_row, Bt_row) 
        
    } # end of Gibbs sampler
    
    
    
    # discard the first n_burn samples
    u_samps = data.frame(u_samps[-c(1:n_burn),]) 
    # u_samps = data.frame(u_samps)
    
    
    return(list(Y = Y, X = X, XtX = XtX, Xty = Xty,
                gibbs_mean = c(Amu, Bt), u_samps = u_samps))
    
} # end of sampleRRR() function ------------------------------------------------



p = 8             # number of columns in X
q = 6             # number of columns in Y
r = 2             # number of columns in B and A
D = r * p + q * r # dimension of each MCMC sample
n = 100           # number of rows in X and Y
sig2 = 1          # fixed for now.
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













