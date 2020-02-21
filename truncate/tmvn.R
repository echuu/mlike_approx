
install.packages("TruncatedNormal")
library(TruncatedNormal)
library(mvtnorm)
library(dplyr)

# path for lenovo
LEN_PATH  = "C:/Users/ericc/mlike_approx"
setwd(LEN_PATH)

# path for dell
# DELL_PATH = "C:/Users/chuu/mlike_approx"
# setwd(DELL_PATH)

# files must be loaded in this order, since *_helper.R file will sometimes
# overwrite a functino defined in hybrid_approx.R depending on the example

source("partition/partition.R")         # load partition extraction functions
source("hybrid_approx.R")               # load main algorithm functions
source("truncate/tmvn_helper.R")        # load psi(), lambda() function
source('extractPartition.R')            # load extractPartition() function

set.seed(1)

N         = 1000
D         = 100

lb        = rep(0, D)
ub        = rep(Inf, D)

mu        =  rep(0, D)    
rho       =  0.2
a         =  1 / (1 - rho)
b         =  - rho * a / ((1 - rho) + D * rho)
omega     =  (1 - rho) * diag(1, D) + rho * outer(rep(1, D), rep(1, D))
omega_inv =  a * diag(1, D) + b * outer(rep(1, D), rep(1, D))
ld_omega  =  ((D - 1) * log(1 - rho) + log((1 - rho) + D * rho)) # logdet(omega)


prior = list(Sigma = omega, Sigma_inv = omega_inv, ld_sigma = ld_omega, 
             rho = rho, D = D)

N_approx = 1
J = N
u_samp = rtmvnorm(N, mu, omega, lb, ub) %>% data.frame


u_df = preprocess(samples, D, prior)

hml_approx = hml(N_approx, D, u_df, J, prior)
hml_approx$hybrid_vec

u_df %>% head
plot(u_samp)





# hml_approx$const_vec    
# hml_approx$taylor_vec
# fun <- function(x, y) (2*pi)^(-D/2) * 1/sqrt(det(Sigma)) * exp(-1/2*(x^2+y^2))
# library(pracma)
# result = integral2(fun, 0, 15, 0, 15, reltol = 1e-50)
# log(result$Q) # -1.223014 for n = 1000
# log(1/2^D)















