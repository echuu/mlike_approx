


library(mvtnorm)
library(dplyr)


sampleRRR = function(n_samps, n_burn, A_0, B_0, p, q, r_0, r, D, N, sig2, del,
                     RANK1 = FALSE) {
    
    # identity matrices
    I_q = diag(1, q)
    I_r = diag(1, r)
    
    eps = matrix(rnorm(N * q, 0, sqrt(sig2)), N, q)   # (n x q) error
    
    
    if (RANK1) {
        X = diag(1, p) # p = N, take X to be the (p x p) identity matrix
    } else {
        X = rmvnorm(N, mean = rep(0, p), diag(1, p))      # (n x p) response
    }
    
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















