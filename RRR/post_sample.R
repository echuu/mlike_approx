
library(mvtnorm)

## likelihood parameters

# ------------------------------------------------------------------------------

# TODO: move all this outside out of the function when finished, should belong 
# ouside of simulation loops


p = 8             # number of columns in X
q = 6             # number of columns in Y
r = 2             # number of columns in B and A
D = r * p + q * r # dimension of each MCMC sample
n = 100           # number of rows in X and Y
sig2 = 10^(-2)    # fixed for now.
del = 10^(-2)     # prior parameter -- one of these is squared version ?


# identity matrices
I_q = diag(1, q)


## generate and FIX parameters A, B
##     NOTE: -- these will be like beta in the regular regression setting
##     in that it is not randomly generated at the beginning of
##     each iteration. it is fixed for all simulations (over grid of N)


### TODO: verify this sampler is doing the right thing
### in order to the sampler, we can use the same starting points, read them 
### in from matlab, and then run both MCMC chains

set.seed(1)
A_0 = matrix(rnorm(p * r, 0, 1), p, r) # (p x r) matrix
B_0 = matrix(rnorm(q * r, 0, 1), q, r) # (q x r) matrix

# ------------------------------------------------------------------------------

## TODO: read in values of X, eps from MATLAB 
## to compute Y, C, XtX, Xty

## generate covariates -- (n x p) matrix, each row ~ MVN(0, I_p)
X = rmvnorm(n, mean = rep(0, p), diag(1, p))

## generate data
eps = matrix(rnorm(n * q, 0, sqrt(sig2)), n, q)     # (n x q) matrix
Y   = X %*% A_0 %*% t(B_0) + eps                      # (n x q) response matrix
C   = A_0 %*% t(B_0)                                # (p x q)

XtX = t(X) %*% X
Xty = t(X) %*% Y


## all things ABOVE should be read in from MATLAB sampler
## we don't include the starting points b/c hopefully the gibbs sampler isn't
## that sensitive to start point. 


# ------------------------------------------------------------------------------

## TODO:
## using matlab RNG parameters above, run the MCMC chain and compare the samples
## at the end of nMCMC runs of the gibbs sampler


## gibbs sampling to obtain samples from (B^T, A) | y, X -----------------------
set.seed(1)
nMCMC = 1500
B = matrix(rnorm(q * r, 0, 1), q, r) # (q x r) matrix for starting point of MCMC
A = matrix(rnorm(p * r, 0, 1), p, r) # (p x r) matrix for starting point of MCMC


u_df = matrix(0, nMCMC, D)
for (g in 1:nMCMC) {
    
    # sample from the conditional distribution: B' | - 
    Btvarpart = sig2 * solve(t(A) %*% XtX %*% A + del^2)
    Btvar = I_q %x% Btvarpart
    Btmu = Btvarpart %*% t(A) %*% Xty / sig2
    
    Bt_row = c(rmvnorm(1, c(Btmu), Btvar)) # (1 x rq) row vector
    
    Mt = solve(t(B) %*% B) %*% (t(B) %*% t(Y) %*% X) %*% solve(XtX)
    M = t(Mt)
    
    BtB = t(B) %*% B
    
    # sample from the conditional distribution A | -
    Avarpart = 1 / sig2 * (BtB %x% XtX)
    Avar = solve(del^2 / sig2 * diag(1, nrow(Avarpart)) + Avarpart)
    Amu = Avar %*% Avarpart %*% c(M)
    
    A_row = c(rmvnorm(1, Amu, Avar))
    
    u_df[g,] = c(A_row, Bt_row)
    
}


# recover matrix A
(A_g = matrix(u_df[g, 1:(p * r)], p, r))

# recover matrix Bt
(Bt_g = matrix(u_df[g,(p * r + 1):D], r, q))


# ------------------------------------------------------------------------------

