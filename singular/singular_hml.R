

library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio


# use stan to draw from the posterior distribution -----------------------------
# setwd("C:/Users/ericc/mlike_approx/singular")
setwd("C:/Users/chuu/mlike_approx/singular")
source("singular_helper.R")
source("C:/Users/chuu/mlike_approx/partition/partition.R")

J         = 2000         # number of MC samples per approximation
N_approx  = 100          # number of approximations
burn_in   = 1000         # number of burn in draws, must be > (N_approx * J)
n_chains  = 4            # number of markov chains to run
stan_seed = 123          # seed

# 6000 iterations needed if we want 2000 MC samples for each approximation
J_iter = 1 / n_chains * N_approx * J + burn_in 


N = 10000  # (pseudo) sample size
D = 2   # dimension of parameter

gamma_dat = list(N = N)

# should give us J * N_approx draws
gamma_fit = stan(file    = 'gamma_sample.stan', 
                 data    = gamma_dat,
                 iter    = J_iter,
                 warmup  = burn_in,
                 chains  = n_chains,
                 seed    = stan_seed,
                 control = list(adapt_delta = 0.99),
                 refresh = 0) 

# to test functionality of preprocess(), don't run the code under this ---------
# manually extract just to see if we're sampling the right thing
u_samp = rstan::extract(gamma_fit, pars = c("u"), permuted = TRUE)

# line below can be called directly after stan() functin call
u_df_all = preprocess(gamma_fit, D, N)

u_df_all %>% head

plot(u_df_all[,1:2])


# u_post = u_samp$u %>% data.frame() # (J * N_approx) x 2

# plot(u_post)

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

    print(paste('iter = ', i, ' -- sampling data for N = ', N, sep = ''))
    
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

## GO TO "START HERE" SECTION --------------------------------------------------


# plot 2-d distributions for a given N
N = 50
filename = paste("data_big/u_df_", N, ".csv", sep = '')
u_df_N = read.csv(filename)
plot(u_df_N[,1:2], main = "N = 50")

# contour plot of 2-d distribution for a given N
ggplot(u_df_N, aes(u1, u2)) + geom_density_2d()


# ------------------------------------------------------------------------------

## test for single N

# read in u_df data
N = 50
filename = paste("u_df_", N, ".csv", sep = '')
u_df_N = read.csv(paste("data/", filename, sep = ''))

# stan_approx will N_approx - dim vector of approximations of the LIL
stan_approx = approx_lil_stan(N_approx, D, N, u_df_N, J)

stan_approx

mean(stan_approx, na.rm = TRUE) # 
var(stan_approx, na.rm = TRUE)  # 

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

approx_mean = colMeans(approx_N)

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




