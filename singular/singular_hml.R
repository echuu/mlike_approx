

library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio


# use stan to draw from the posterior distribution -----------------------------
setwd("C:/Users/ericc/mlike_approx/singular")
source("singular_helper.R")

J         = 2000         # number of MC samples per approximation
N_approx  = 10           # number of approximations
burn_in   = 1000         # number of burn in draws, must be > (N_approx * J)
n_chains  = 4            # number of markov chains to run
stan_seed = 123          # seed

J_iter = 1 / n_chains * N_approx * J + burn_in 


N = 50  # (pseudo) sample size
D = 2   # dimension of parameter

gamma_dat = list(N = N)

# should give us J * N_approx draws
gamma_fit = stan(file    = 'gamma_sample.stan', 
                 data    = gamma_dat,
                 iter    = J_iter,
                 warmup  = burn_in,
                 chains  = n_chains,
                 seed    = stan_seed,
                 control = list(adapt_delta = 0.99)) 

# to test functionality of preprocess(), don't run the code under this ---------
# manually extract just to see if we're sampling the right thing
u_samp = rstan::extract(gamma_fit, pars = c("u"), permuted = TRUE)

u_post = u_samp$u %>% data.frame() # (J * N_approx) x 2

plot(u_post)

# ------------------------------------------------------------------------------

# line below can be called directly after stan() functin call
u_df_all = preprocess(gamma_fit, D, N)

plot(u_df_all[,1:2])


filename = paste("u_df_", N, ".csv", sep = '')
test_write = write.csv(u_df_all, filename, row.names = FALSE)
test_read = read.csv("u_df_50.csv") # (J * N_approx) x 3

# ------------------------------------------------------------------------------


# repeat above analysis for a grid of N ----------------------------------------

N_vec  = c(50, 100, 200, 500, 750, 1000, 2000, 4000, 8000, 10000)
N_sims = length(N_vec)

for (i in 1:N_sims) {
    
    
    N = N_vec[i]   # pseudo-sample size
    D = 2          # dimension of parameter
    
    print(paste('iter = ', i, ' -- sampling data for N = ', N, sep = ''))
    
    gamma_dat = list(N = N)
    
    # should give us (J * N_approx) draws
    gamma_fit_N = stan(file    =  'gamma_sample.stan', 
                       data    =  gamma_dat,
                       iter    =  J_iter,
                       warmup  =  burn_in,
                       chains  =  n_chains,
                       seed    =  stan_seed,
                       control =  list(adapt_delta = 0.99)) 
    
    u_df_N = preprocess(gamma_fit_N, D, N)
    
    # save data 
    filename_N = paste("u_df_", N, ".csv", sep = '')
    write.csv(u_df_N, filename_N, row.names = FALSE)

}

u_df_50 = read.csv("u_df_50.csv")

plot(u_df_50[,1:2])


u_df_1000 = read.csv("u_df_10000.csv")

plot(u_df_1000[,1:2])


# ------------------------------------------------------------------------------

ggplot(u_df_all, aes(u1, u2)) + geom_density_2d()


u_df = u_df_all[1:J,]

u_rpart = rpart(psi_u ~ ., u_df)

# (3.1) obtain the (data-defined) support for each of the parameters
param_support = matrix(NA, D, 2) # store the parameter supports row-wise

for (d in 1:D) {
    param_d_min = min(u_df[,d])
    param_d_max = max(u_df[,d])
    
    param_support[d,] = c(param_d_min, param_d_max)
}

# (3.2) obtain the partition --- moment of truth!!
u_partition = paramPartition(u_rpart, param_support)  # partition.R

head(u_partition)

# (3.3) organize all data into single data frame (see partition.R for format)

# extracts u_star, representative point of each partition (u_star \in R^D)
# psi_hat, leaf_id, u1_star, u2_star, ... , uD_star, 
#                   u1_lb, u1_ub, ...uD_lb, uD_ub
param_out = u_star(u_rpart, u_df, u_partition, D)

head(param_out)

# test lambda 
# k = 2
# u_test = c(param_out[k,]$u1_star, param_out[k,]$u2_star)
# exp(-psi(u_test, N))     # (1 x 1) -- c_k[k] calculation
# lambda(u_test, N) # (D x 1) -- lambda(u_star) calculation


## (4) begin main algorithm 
n_partitions = nrow(u_partition)
c_k = numeric(n_partitions)
zhat = numeric(n_partitions)

for (k in 1:n_partitions) {
    
    # u_star_k = (mu_k, sigma_sq_k)
    # c_k[k] = exp(-psi(param_out[k,]$u1_star, 
    #                   param_out[k,]$u2_star,
    #                   y, m_0, w_0, r_0, s_0)) # (1 x 1)
    
    # TODO: avoid explicit definition of each representative point
    u = c(param_out[k,]$u1_star, param_out[k,]$u2_star)
    
    c_k[k] = exp(-psi(u, N)) # (1 x 1)
    
    l_k = lambda(u, N)
    
    integral_d = numeric(D) # store each component of the D-dim integral 
    
    # nothing to refactor in this loop (i think?) since we're just iterating
    # thru each of the integrals and computing an exponential term
    for (d in 1:D) {
        
        # verify these -- these need to be recalculated if the form of param_out
        # changes (if columns get shuffled)
        
        # col id will change for D > 2
        # DONE: generalize this better so there's less obscure calculation
        # col_id_lb = 5 + 2 * (d - 1)
        # col_id_ub = col_id_lb + 1
        
        # updated 1/14: find column id of the first lower bound
        col_id_lb = grep("u1_lb", names(param_out))
        col_id_ub = col_id_lb + 1
        
        # d-th integral computed in closed form
        integral_d[d] = - 1 / l_k[d] * 
            exp(- l_k[d] * (param_out[k, col_id_ub] - param_out[k, col_id_lb]))        
        
    }
    
    zhat[k] = prod(c_k[k], integral_d)
}

log(sum(zhat))




