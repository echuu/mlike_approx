
##  RRR_sample.R 
##  read in the posterior samples written to .csv from MATLAB script

library(dplyr)

# path for lenovo
LEN_PATH  = "C:/Users/ericc/mlike_approx"
setwd(LEN_PATH)

# files must be loaded in this order, since *_helper.R file will sometimes
# overwrite a functino defined in hybrid_approx.R depending on the example

source("partition/partition.R")         # load partition extraction functions
source("hybrid_approx.R")               # load main algorithm functions
source("mvn/mvn_helper.R")              # load psi(), lambda() function
source('extractPartition.R')            # load extractPartition() function


## read in true parameters generated/computed via MATLAB script ----------------
##
# covariates
# read in design matrix
X = read.csv("RRR/X.csv", header = FALSE) %>% as.matrix

# read in parameters A, B
A = read.csv("RRR/A.csv", header = FALSE) %>% as.matrix
B = read.csv("RRR/B.csv", header = FALSE) %>% as.matrix

# read in error: eps
eps = read.csv("RRR/eps.csv", header = FALSE) %>% as.matrix

# compute: Y = X * A * B' + eps
Y = X %*% A %*% t(B) + eps

# compute C = A * B'
C = A %*% t(B)

##
## -----------------------------------------------------------------------------

# dimensions can be computed from the dimensions of the input above

# various dimension settings
p = ncol(X)    # number of columns in X
q = ncol(Y)    # number of columns in Y
r = ncol(A)    # number of columns in B and A
n = nrow(X)    # number of rows in X and Y


# prior parameters
sig2 = 10^(-2);  # fixed for now
del  = 10^(-2);

# dimension of the parameter, u (drawn from posterior)
d = r * p + r * q

# store parameters in param object
param_list = list(p = p, q = q, r = r, n = n, d = d,   # dimensions variables
                  Y = Y, X = X,                        # response, design matrix
                  sig2 = sig2, del = del)              # prior params


burn = 500   # number of samples to discard
J = 200      # number of MCMC samples to use per approximation

# read in posterior samples 
u_samps = read.csv("RRR/u_df_rrr.csv", header = FALSE)
dim(u_samps)

### can stick the following directly into preprocess function
u_samps = u_samps[-c(1:burn),]
u_samps_sub = u_samps[1:J,]    # (J x d) -- MCMC samples stored row-wise


## evaluate psi() function for each of the posterior samples (row-wise)
u_df = preprocess(u_samps_sub, d, param_list)


# for now, when moving this into helper function, just load the helper function 
# after loading hybrid_approx.R so that it overrides it
# this function, and psi() both use 'params' arg instead of 'prior' --
# should be switching over to params anyway, it's not always just the prior 
# parameters getting passed into preprocess, psi..
preprocess = function(post_samps, D, params) {
    
    psi_u = apply(post_samps, 1, psi, params = params) %>% unname() # (J x 1)
    
    # (1.2) name columns so that values can be extracted by partition.R
    u_df_names = character(D + 1)
    for (d in 1:D) {
        u_df_names[d] = paste("u", d, sep = '')
    }
    u_df_names[D + 1] = "psi_u"
    
    # populate u_df
    u_df = cbind(post_samps, psi_u) # J x (D + 1)
    names(u_df) = u_df_names
    
    
    return(u_df)
    
} # end of preprocess() function


psi = function(u, params) {
    
    # note that the norm(,type = 'F') runs 10x faster than computing the trace
    
    n = params$n # num of rows in X, Y
    p = params$p # num rows of A
    r = params$r # num cols of A, num rows of B^T
    q = params$q # num cols of B^T
    
    sig2 = params$sig2
    del  = params$del
    
    Y = params$Y
    X = params$X
    
    # extract and reshape u to get the matrices: A, B^T
    A_post  = matrix(u[1:(p * r)], p, r)
    Bt_post = matrix(u[(p * r + 1):d], r, q)
    
    const_term = - 0.5 * n * q * log(2 * pi * sig2)
    exp_term   = -1/(2 * sig2) * 
        norm(Y - X %*% A_post %*% Bt_post, type = 'F')^2 + 
        del^2 / (2 * sig2) * 
        (norm(Bt_post, type = 'F')^2 + norm(A_post, type = 'F')^2)
    
    return(const_term + exp_term)
    
} # end of psi() function
    

lambda = function(u, params) {
    
    
    
    
} # end of lambda() function
 






# library(microbenchmark)
# microbenchmark("frob" = { a = norm(B, type = 'F')^2 },
#                "trace" = { b = sum(diag((t(B) %*% B))) })


## evaluate lambda() function

































## uncomment the R code below to check that the values match matlab output
# obtain A, B^T
# A_g  = reshape(u_df(g,1:p*r), p, r);    # recover the matrix A
# Bt_g = reshape(u_df(g,p*r+1:D), r, q);  # recover the matrix B^T =: B

# # A is (p x r)
# matrix(u_df[1000, 1:(p*r)], p, r) # fill the matrix COLUMN-wise (what we want)
# 
# # B^T is (r x q)
# matrix(u_df[1000, (p*r+1):d], r, q)
