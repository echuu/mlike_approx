

library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio


# use stan to draw from the posterior distribution -----------------------------
# setwd("C:/Users/ericc/mlike_approx/singular")
setwd("C:/Users/ericc/mlike_approx/singular")
source("singular_helper.R")
source("C:/Users/ericc/mlike_approx/partition/partition.R")

J         = 2000         # number of MC samples per approximation
N_approx  = 10           # number of approximations
burn_in   = 1000         # number of burn in draws, must be > (N_approx * J)
n_chains  = 4            # number of markov chains to run
stan_seed = 123          # seed

# 6000 iterations needed if we want 2000 MC samples for each approximation
J_iter = 1 / n_chains * N_approx * J + burn_in 

# values of N for which we will compute + approximate the LIL
N_vec = seq(50, 5000, 10)

N_vec = c(50, 100)

N_sims = length(N_vec)

for (i in 1:N_sims) {
    
    
    N = N_vec[i]   # pseudo-sample size
    
    print(paste('iter = ', i, ' -- sampling data for N = ', N, sep = ''))
    
    # repeat the below code K times
    gamma_dat = list(N = N)
    
    # should give us (J * N_approx) draws
    gamma_fit_N = stan(file    =  'gamma_sample.stan', 
                       data    =  gamma_dat,
                       iter    =  J_iter,
                       warmup  =  burn_in,
                       chains  =  n_chains,                  # num chains used
                       seed    =  stan_seed,
                       control =  list(adapt_delta = 0.99),  # small step size
                       refresh = 0)                          # suppress output
    
    u_df_N = preprocess(gamma_fit_N, D, N)
    
    
    # save data 
    filename_N = paste("u_df_", N, ".csv", sep = '')
    write.csv(u_df_N, filename_N, row.names = FALSE)
    print(paste('iter = ', i, ' -- writing out ', filename_N, sep = ''))
    
}



## START HERE ------------------------------------------------------------------
## get approximation for vector of N \in N_vec
## store approximations for each N col-wise in (N_approx) x (N_sims) matrix
approx_N = matrix(NA, N_approx, N_sims) # (N_approx x n_cases)

for (i in 1:N_sims) {
    
    N = N_vec[i]   # pseudo-sample size
    
    filename = paste("u_df_", N, ".csv", sep = '')
    u_df_N = read.csv(paste("data_big/", filename, sep = ''))
    
    print(paste('iter = ', i, ' -- calculating LIL for N = ', N, " (", 
                N_approx, " approximations)", sep = ''))
    
    approx_N[,i] = approx_lil_stan(N_approx, D, N, u_df_N, J)
}


approx_med = apply(approx_N, 2, median) %>%  unlist() %>% unname() 
approx_med

approx_sds = apply(approx_N, 2, sd) %>%  unlist() %>% unname() 
approx_sds

# save simulation results
write.csv(approx_N, "singular_lil.csv", row.names = F)

test_read = read.csv("singular_lil.csv")

# asymptotic analysis ----------------------------------------------------------

z_n = exp(approx_mean)
log_z_n = approx_mean
log_n = log(N_vec)

lil_df = data.frame(z_n, log_z_n, log_n)

ggplot(lil_df, aes(log_n, log_z_n)) + geom_point()

lm(log_z_n ~ log_n, lil_df)











