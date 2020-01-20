

library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio


# use stan to draw from the posterior distribution -----------------------------
# setwd("C:/Users/ericc/mlike_approx/singular")
source("C:/Users/ericc/mlike_approx/partition/partition.R")
setwd("C:/Users/ericc/mlike_approx/singular")
source("singular_helper.R")



# STAN SETTINGS ----------------------------------------------------------------
J         = 2000         # number of MC samples per approximation
N_approx  = 50           # number of approximations
burn_in   = 2000         # number of burn in draws
n_chains  = 4            # number of markov chains to run
stan_seed = 123          # seed

J_iter = 1 / n_chains * N_approx * J + burn_in 
# ------------------------------------------------------------------------------


# values of N for which we will compute + approximate the LIL
N_vec = seq(50, 5000, 10)
N_vec = c(22026)   # pseudo sample size


N_vec_log = seq(4, 10, 0.05)        # sample size that is uniform over log scale
N_vec     = floor(exp(N_vec_log))   # sample size to use to generate data
D = 2                               # dimension of parameter


# store approximations corresponding to each sample size
approx_N = matrix(NA, N_approx, length(N_vec))

for (i in 1:length(N_vec)) {
    
    
    N = N_vec[i]   # pseudo-sample size
    
    # print(paste('iter = ', i, ' -- sampling data for N = ', N, sep = ''))
    
    gamma_dat = list(N = N)
    
    # should give us (J * N_approx) draws
    gamma_fit_N = stan(file    =  'gamma_sample.stan', 
                       data    =  gamma_dat,
                       iter    =  J_iter,
                       warmup  =  burn_in,
                       chains  =  n_chains,                  
                       seed    =  stan_seed,
                       control =  list(adapt_delta = 0.99),  
                       refresh = 0)                          
    
    u_df_N = preprocess(gamma_fit_N, D, N)
    
    print(paste('iter = ', i, ' -- calculating LIL for N = ', N, 
                ' (', N_approx, ' approximations)', sep = ''))
    
    
    # N_approx, J settings indicate that J MC samples will use in each of the
    # N_approx estiamtes of the LIL -> return vector of length N_approx
    approx_N[,i] = approx_lil_stan(N_approx, D, N, u_df_N, J)
    
}

approx_N

write.csv(approx_N, "approx_N50_J2000_grid121.csv", row.names = F)
test_read = read.csv("approx_N50_J2000_grid121.csv")


colMeans(approx_N)











