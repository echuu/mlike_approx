

library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio
library("ggpmisc")

# use stan to draw from the posterior distribution -----------------------------

# path for lenovo
# LEN_PATH  = "C:/Users/ericc/mlike_approx"
# setwd(LEN_PATH)

# path for dell
DELL_PATH = "C:/Users/chuu/mlike_approx"
setwd(DELL_PATH)

source("partition/partition.R")         # load partition extraction functions
source("hybrid_approx.R")               # load main algorithm functions

source("singular/singular_helper.R")    # load psi(), lambda()



# STAN SETTINGS ----------------------------------------------------------------
J         = 1000         # number of MC samples per approximation
N_approx  = 10           # number of approximations
burn_in   = 2000         # number of burn in draws
n_chains  = 4            # number of markov chains to run
stan_seed = 123          # seed

J_iter = 1 / n_chains * N_approx * J + burn_in 
# ------------------------------------------------------------------------------

# GLOBAL MODEL SETTINGS --------------------------------------------------------
D = 2                    # dimension of parameter

# one run of the algorithm -----------------------------------------------------
set.seed(1)
N = 1000

gamma_dat = list(N = N) # for STAN sampler
prior     = list(N = N) # for evaluation of psi, lambda

# (1) generate posterior samples -- should give us (J * N_approx) draws
gamma_fit_N = stan(file    =  'singular/gamma_sample.stan', 
                   data    =  gamma_dat,
                   iter    =  J_iter,
                   warmup  =  burn_in,
                   chains  =  n_chains,                  
                   seed    =  stan_seed,
                   control =  list(adapt_delta = 0.99),  
                   refresh = 0)           

u_samp = rstan::extract(gamma_fit_N, pars = c("u"), permuted = TRUE)
u_post = u_samp$u %>% data.frame() # (J * N_approx) x 2

# (2) evaluate posterior samples using psi(u)
u_df_N = preprocess(u_post, D, prior)

ggplot(u_df_N, aes(u1, u2)) + geom_point()

# (3) run algorithm to obtain N_approx estimates of the LIL
approx = approx_lil(N_approx, D, u_df_N, J, prior)

approx

mean(approx) # 9.044141 for D = 2, N = 1000

# ------------------------------------------------------------------------------


# run algorithm over grid of N -------------------------------------------------


# values of N for which we will compute + approximate the LIL

N_vec_log = seq(4, 10, 0.05)       # sample size that is uniform over log scale
N_vec     = floor(exp(N_vec_log))  # sample size to use to generate data


N_vec     = c(1000)                # sample size to use to generate data
set.seed(1)

# store approximations corresponding to each sample size
approx_N = matrix(NA, N_approx, length(N_vec))

for (i in 1:length(N_vec)) {
    
    
    N = N_vec[i]   # pseudo-sample size
    
    # print(paste('iter = ', i, ' -- sampling data for N = ', N, sep = ''))
    
    gamma_dat = list(N = N) # for STAN sampler
    prior     = list(N = N) # for evaluation of psi, lambda
    
    # (1) generate posterior samples -- should give us (J * N_approx) draws
    gamma_fit_N = stan(file    =  'singular/gamma_sample.stan', 
                       data    =  gamma_dat,
                       iter    =  J_iter,
                       warmup  =  burn_in,
                       chains  =  n_chains,                  
                       seed    =  stan_seed,
                       control =  list(adapt_delta = 0.99),  
                       refresh = 0)           
    
    u_samp = rstan::extract(gamma_fit_N, pars = c("u"), permuted = TRUE)
    u_post = u_samp$u %>% data.frame() # (J * N_approx) x 2
    
    # (2) evaluate posterior samples using psi(u)
    u_df_N = preprocess(u_post, D, prior)
    
    # (3) run algorithm to obtain N_approx estimates of the LIL
    
    print(paste('iter = ', i, ' -- calculating LIL for N = ', N, 
                ' (', N_approx, ' approximations)', sep = ''))
    
    
    # N_approx, J settings indicate that J MC samples will use in each of the
    # N_approx estiamtes of the LIL -> return vector of length N_approx
    approx_N[,i] = approx_lil(N_approx, D, u_df_N, J, prior)
    
}

approx_N

colMeans(approx_N)

# write.csv(approx_N, "approx_N50_J2000_grid121.csv", row.names = F)
# test_read = read.csv("approx_N50_J2000_grid121.csv")


log_z_n = colMeans(approx_N)
log_n   = log(N_vec)

lil_df = data.frame(log_z_n, log_n)

formula1 = y ~ x

ggplot(lil_df, aes(log_n, log_z_n)) + geom_point() + 
    labs(title = "50 approximations (1000 MC samples each) for each N") + 
    geom_smooth(method = lm, se = T, formula = formula1) +
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16)

lm(log_z_n ~ log_n, lil_df) # slope should be -0.25











