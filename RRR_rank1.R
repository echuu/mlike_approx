
## RRR_rank1.R -----------------------------------------------------------------

library(dplyr)
library(mvtnorm)


# path for lenovo
LEN_PATH  = "C:/Users/ericc/mlike_approx"
setwd(LEN_PATH)

# files must be loaded in this order, since *_helper.R file will sometimes
# overwrite a function defined in hybrid_approx.R depending on the example

source("partition/partition.R")         # load partition extraction functions
source("hybrid_approx_v1.R")            # load main algorithm functions
# source("RRR/RRR_rank1_helper.R")        # load psi(), lambda(), preprocess()
source("RRR/RRR_C_helper.R")            # load psi(), lambda(), preprocess()
source("RRR/gibbs_RRR.R")               # load sampleRRR()
source('extractPartition.R')            # load extractPartition() function


# ------------------------------------------------------------------------------

q = 3
p = 2
r = 1
D = r * p + r * q
D_C = p * q        # dimension of the integral
N = 100

sig2 = 1           # fixed for now.
del = 10^(-2)     # prior parameter -- one of these is squared version ?

## model : Y = ab' + E
# Y is (n x q)
# X is (n x p)
# a is (p x 1)
# b is (q x 1)
# E is (p x q)

# ------------------------------------------------------------------------------
set.seed(1)
A_0 = matrix(rnorm(p * r, 0, 1), p, r) # (p x r) matrix
B_0 = matrix(rnorm(q * r, 0, 1), q, r) # (q x r) matrix
C_0 = A_0 %*% t(B_0)


nMCMC = 300       # number of MCMC samples from the posterior AFTER burnin
nBurn = 500       # number of samples to discard

gibbs_obj = sampleRRR(nMCMC, nBurn, A_0, B_0, p, q, r, r, D, N, sig2, del)

param_list = list(p = p, q = q, r = r, n = N, d = D,          # dimensions vars
                  Y = gibbs_obj$Y, X = gibbs_obj$X,           # response, design
                  XtX = gibbs_obj$XtX, Xty = gibbs_obj$Xty,   # precompute
                  sig2 = sig2, del = del)                     # prior params

u_samps = gibbs_obj$u_samps
c_samps = gibbs_obj$c_samps

u = u_samps[1,]
a_post = u[1:p] %>% unlist %>% unname # (p x 1)
b_post = u[-(1:p)] %>% unlist %>% unname

lambda(unlist(unname(c_samps[1,])), param_list)

a_post %*% t(b_post)

matrix(c_samps[1,], p, q)

head(u_samps)
# ------------------------------------------------------------------------------

# u_df = preprocess(u_samps, D, param_list) # nMCMC x (d + 1) 

c_df = preprocess(c_samps, D_C, param_list)

hml_approx = hml(1, D_C, c_df, nMCMC, param_list)
hml_approx$hybrid_vec
hml_approx$n_taylor
hml_approx$n_const
hml_approx$taylor_vec


library(matrixcalc)

marglik = function(u) {
    
    a = u[1:p]
    b = tail(u, q)
    
    # (-N * q / 2) * log(2 * pi * sig2) -
    #     1 / (2 * sig2) * 
    #     norm(gibbs_obj$Y - gibbs_obj$X %*% a %*% t(b), type = 'F')^2
    
    (2 * pi * sig2)^(-N * q / 2) *
         exp(-1/(2 * sig2) *
                 norm(gibbs_obj$Y - gibbs_obj$X %*% a %*% t(b), type = 'F')^2)

    # (2 * pi * sig2)^(-p * q/2) *
    #     exp(-1/(2 * sig2) * norm(rank1_obj$Y - a %*% t(b), type = 'F')^2) *
    #     (2 * pi * sig2)^(-(p + q) / 2) * (del^2)^((p + q) / 2) *
    #     exp(-del^2 / (2 * sig2) * (t(a) %*% a + t(b) %*% b))
    
}

stable_marglik = function(u) {
    
    a = u[1:p]
    b = tail(u, q)
    
    c_mat = a %*% t(b)
    
    exp(-1/(2 * sig2) * 
            matrix.trace(t(c_mat) %*% gibbs_obj$XtX %*% c_mat - 
                             2 * t(gibbs_obj$Xty) %*% c_mat))
    
}

library(cubature)
marglik_numer = adaptIntegrate(marglik, 
                               lowerLimit = c(-Inf, -Inf, -Inf, -Inf, -Inf), 
                               upperLimit = c(Inf, Inf, Inf, Inf, Inf),
                               tol = 1e-4)

marglik_numer$integral
log(marglik_numer$integral)  # -423.2818

stable_numer = adaptIntegrate(stable_marglik, 
                              lowerLimit = c(-Inf, -Inf, -Inf, -Inf, -Inf), 
                              upperLimit = c(Inf, Inf, Inf, Inf, Inf),
                              tol = 1e-4)


(scaled_marglik = -N * q / 2 * log(2 * pi * sig2) - 
    1 / (2 * sig2) * matrix.trace(t(gibbs_obj$Y) %*% gibbs_obj$Y) + 
    log(stable_numer$integral))



# testing for preventing underflow ---------------------------------------------


u = u_samps[1,]
a_post = u[1:p] %>% unlist %>% unname # (p x 1)
b_post = u[-(1:p)] %>% unlist %>% unname
c_post = a_post %*% t(b_post)

# these quantities should be equal ---------------------------------------------
log((2 * pi * sig2)^(-N * q / 2) *
    exp(-1/(2 * sig2) *
            norm(gibbs_obj$Y - gibbs_obj$X %*% a_post %*% t(b_post), type = 'F')^2))

-N * q / 2 * log(2 * pi * sig2) - 
    1 / (2 * sig2) * matrix.trace(t(gibbs_obj$Y) %*% gibbs_obj$Y) -
    1 / (2 * sig2) * matrix.trace(t(c_post) %*% gibbs_obj$XtX %*% c_post - 2 * t(gibbs_obj$Xty) %*% c_post)
    





# ------------------------------------------------------------------------------


a_0 = rnorm(p) # (p x 1)
b_0 = rnorm(q) # (q x 1)
c_0 = a_0 %*% t(b_0)
lil_hat = numeric(50)

for(i in 1:50) {
    
    nMCMC = 500       # number of MCMC samples from the posterior AFTER burnin
    nBurn = 1000       # number of samples to discard
    
    rank1_obj = sampleRRR(nMCMC, nBurn, a_0, b_0, p, q, r, r, D, N, sig2, del, T)
    
    u_samps = rank1_obj$u_samps # (nMCMC x D)
    
    param_list = list(p = p, q = q, r = r, n = N, d = D,  # dimensions vars
                      Y = rank1_obj$Y, X = rank1_obj$X,   # response, design
                      XtX = rank1_obj$XtX, Xty = rank1_obj$Xty,
                      sig2 = sig2, del = del)             # prior params
    
    # evaluate psi(u) for each of the posterior samples
    u_df = preprocess(u_samps, D, param_list) # nMCMC x (d + 1) 
    
    # generate hybrid approximation
    hml_approx = hml(1, D, u_df, nMCMC, param_list)
    lil_hat[i] = hml_approx$hybrid_vec
    
    
}

mean(lil_hat)
sd(lil_hat)


a_0 = rnorm(p) # (p x 1)
b_0 = rnorm(q) # (q x 1)

c_0 = a_0 %*% t(b_0)

nMCMC = 500       # number of MCMC samples from the posterior AFTER burnin
nBurn = 1000       # number of samples to discard

rank1_obj = sampleRRR(nMCMC, nBurn, a_0, b_0, p, q, r, r, D, N, sig2, del, T)

u_samps = rank1_obj$u_samps # (nMCMC x D)

u = u_samps[300,] %>% unname %>% unlist

a_post = u[1:p]      # (p x 1)
b_post = tail(u, q)  # (q x 1)

(c_post = a_post %*% t(b_post))

c_0
# ------------------------------------------------------------------------------


# param_list = list(p = p, q = q, r = r, n = N, d = D,  # dimensions vars
#                   Y = rank1_obj$Y, X = rank1_obj$X,   # response, design
#                   XtX = rank1_obj$XtX, Xty = rank1_obj$Xty,
#                   sig2 = sig2, del = del)             # prior params
# 
# # evaluate psi(u) for each of the posterior samples
# u_df = preprocess(u_samps, D, param_list) # nMCMC x (d + 1) 
# 
# # generate hybrid approximation
# hml_approx = hml(1, D, u_df, nMCMC, param_list)
# hml_approx$hybrid_vec
# hml_approx$n_taylor


# compute the numerical integral
library(pracma)
fun <- function(x, y) exp(-n*x^2*y^4)
result = integral2(fun, 0, 1, 0, 1, reltol = 1e-50)
log(result$Q) # -1.223014 for n = 1000



library(cubature)

f = function(x) {
    a = x[1:2]
    b = x[3:4]
    2/3 * (a[1] * b[1] + a[2] * b[2]) 
}

f = function(x) {
    a = x[1:2]
    b = x[3:4]
    2/3 * (t(a) %*% b) 
}

marglik = function(u) {
    
    a = u[1:p]
    b = tail(u, q)
    
    (2 * pi * sig2)^(-p * q/2) * 
        exp(-1/(2 * sig2) * norm(rank1_obj$Y - a %*% t(b), type = 'F')^2) *
        (2 * pi * sig2)^(-(p + q) / 2) * (del^2)^((p + q) / 2) * 
        exp(-del^2 / (2 * sig2) * (t(a) %*% a + t(b) %*% b))
    
}


adaptIntegrate(marglik, 
               lowerLimit = c(-Inf, -Inf, -Inf, -Inf, -Inf), 
               upperLimit = c(Inf, Inf, Inf, Inf, Inf))


marglik_numer = adaptIntegrate(marglik, 
                               lowerLimit = c(-Inf, -Inf, -Inf, -Inf, -Inf), 
                               upperLimit = c(Inf, Inf, Inf, Inf, Inf))

log(marglik_numer$integral)

r1_logML = 0.00000169861
log(r1_logML)
















