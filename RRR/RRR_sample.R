
##  RRR_sample.R 
##  read in the posterior samples written to .csv from MATLAB script

library(dplyr)
library(mvtnorm)


# path for lenovo
LEN_PATH  = "C:/Users/ericc/mlike_approx"
setwd(LEN_PATH)

# files must be loaded in this order, since *_helper.R file will sometimes
# overwrite a functino defined in hybrid_approx.R depending on the example

source("partition/partition.R")         # load partition extraction functions
source("hybrid_approx.R")               # load main algorithm functions
source("RRR/RRR_helper.R")              # load psi(), lambda(), preprocess()
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


# various dimension settings
# dimensions can be computed from the dimensions of the input above
p = ncol(X)    # number of columns in X
q = ncol(Y)    # number of columns in Y
r = ncol(A)    # number of columns in B and A
n = nrow(X)    # number of rows in X and Y

# prior parameters
sig2 = 10^(-2);  # fixed for now
del  = 10^(-2);

# dimension of the parameter, u (drawn from posterior)
d = r * p + r * q

# various identity matrices
I_p = diag(1, p)
I_q = diag(1, q)
I_n = diag(1, n)
I_r = diag(1, r)

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
# u_samps = u_samps[-c(1:burn),]
# u_samps_sub = u_samps[1:J,]    # (J x d) -- MCMC samples stored row-wise


## evaluate psi() function for each of the posterior samples (row-wise)
u_df = preprocess(u_samps, d, param_list) # J x (d + 1) 

source("hybrid_approx_v1.R")               # load main algorithm functions
source("RRR/RRR_helper.R")              # load psi(), lambda(), preprocess()
u_df = preprocess(data.frame(u_df), d, param_list) # J x (d + 1) 

# source("hybrid_approx.R")               # load main algorithm functions

hml_approx = hml(1, d, u_df, J, param_list)
hml_approx$verbose_partition %>% dim
hml_approx$hybrid_vec
hml_approx$ck_3

## -----------------------------------------------------------------------------
u = u_df[1, 1:d] %>% unlist %>% unname
lambda(u, param_list)


# replicate structure of the remaining part of the call to the main algorithm


## -----------------------------------------------------------------------------

J = nrow(u_df) # number of MCMC samples to use in the approximation

source("hybrid_approx_v1.R")                   # load main algorithm functions
hml_approx = hml(1, d, u_df, J, param_list) 


stable_ck3 = hml_approx$ck_3
stable_ck3

old_ck3 = hml_approx$ck_3
old_ck3

cbind(stable_ck3, old_ck3)




## -----------------------------------------------------------------------------

hml_approx$const_vec
hml_approx$taylor_vec
hml_approx$hybrid_vec

hml_approx$n_taylor
hml_approx$verbose_partition
hml_approx$taylor_approx

hml_approx$ck_2

dim(hml_approx$lambda) # 6 partitions x 28 dim parameter

hml_approx$lambda[6,]

hml_approx$ck_3

ind = 20
hml_approx$ck_3[ind]

upper = hml_approx$partition$u20_ub[6]
lower = hml_approx$partition$u20_lb[6]

(l_k_d = hml_approx$lambda[6,ind])

# log(-1 / hml_approx$lambda[6,ind] * 
#     exp(- hml_approx$lambda[6,ind] * hml_approx$partition$u18_ub[6]) - 
#     exp(- hml_approx$lambda[6,ind] * hml_approx$partition$u18_lb[6]))


library(VGAM) # log1mexp(x)


hml_approx$ck_3[ind]

- l_k_d * upper + log(- 1 / l_k_d * (1 - exp(-l_k_d * lower + l_k_d * upper)))

# for lambda > 0
- l_k_d * lower - log(l_k_d) + log1mexp(l_k_d * upper - l_k_d * lower)

# for lambda < 0
-log(-l_k_d) - l_k_d * upper + log1mexp(l_k_d * lower - l_k_d * upper)









# exp(-l_k_d * lower + l_k_d * upper) # overflowing
# exp(-l_k_d * upper + l_k_d * lower)
# 
# 
# exp(-l_k_d * lower) * exp(l_k_d * upper)
# 
# -l_k_d * lower + l_k_d * upper 
# 
# 
# log(1 - exp(-l_k_d * lower + l_k_d * upper))
# 
# 
# library(brms)
# library(VGAM) # log1mexp(x)
# 
# log1p(- exp(-l_k_d * lower + l_k_d * upper))
# 
# # overflowing for exp(x), x > 700
# log1mexp(l_k_d * upper - l_k_d * lower)
# 
# hml_approx$ck_3[ind]








# 
# LIL_const[,k]  = hml_approx$const_vec
# LIL_taylor[,k] = hml_approx$taylor_vec
# LIL_hybrid[,k] = hml_approx$hybrid_vec


# library(microbenchmark)
# microbenchmark("manual" = { a = - 0.5 * N * log(2 * pi * sigmasq) - 
#     1 / (2 * sigmasq) * sum((y - X %*% beta)^2) },
#                "dnorm" = { b = sum(dnorm(y, X %*% mu_beta, sqrt(sigmasq), log = T)) })
# 
# 
# 
# - 0.5 * N * log(2 * pi * sigmasq) - 
#     1 / (2 * sigmasq) * sum((y - X %*% beta)^2)



## -----------------------------------------------------------------------------


# library(microbenchmark)
# microbenchmark("frob" = { a = norm(B, type = 'F')^2 },
#                "trace" = { b = sum(diag((t(B) %*% B))) })



## uncomment the R code below to check that the values match matlab output
# obtain A, B^T
# A_g  = reshape(u_df(g,1:p*r), p, r);    # recover the matrix A
# Bt_g = reshape(u_df(g,p*r+1:D), r, q);  # recover the matrix B^T =: B

# # A is (p x r)
# matrix(u_df[1000, 1:(p*r)], p, r) # fill the matrix COLUMN-wise (what we want)
# 
# # B^T is (r x q)
# matrix(u_df[1000, (p*r+1):d], r, q)
