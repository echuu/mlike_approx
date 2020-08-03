# setup global environment, load in algo functions
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           
source("C:/Users/ericc/mlike_approx/singular/singular_d_helper.R")

library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio
library("ggpmisc")
options(mc.cores = parallel::detectCores()) 

# path for lenovo
stan_sampler = 'C:/Users/ericc/mlike_approx/singular/gamma_d.stan'

# STAN SETTINGS ----------------------------------------------------------------
J         = 2000         # number of MC samples per approximation
N_approx  = 1            # number of approximations
burn_in   = 2000         # number of burn in draws
n_chains  = 4            # number of markov chains to run
stan_seed = 123          # seed

J_iter = 1 / n_chains * N_approx * J + burn_in 

D = 4

# one run of the algorithm -----------------------------------------------------
set.seed(123)

N_vec_log = seq(1, 15, by = 0.1)             # sample size grid unif in log
N_vec     = floor(exp(N_vec_log)) %>% unique # sample size to generate data
print(length(N_vec))               # number of different sample sizes

N = 100

gamma_dat = list(N = N) # for STAN sampler
prior     = list(N = N) # for evaluation of psi, lambda

# (1) generate posterior samples -- should give us (J * N_approx) draws
gamma_fit_N = stan(file    =  stan_sampler, 
                   data    =  gamma_dat,
                   iter    =  J_iter,
                   warmup  =  burn_in,
                   chains  =  n_chains,                  
                   control =  list(adapt_delta = 0.99),  
                   refresh = 0)     

u_samp = rstan::extract(gamma_fit_N, pars = c("u"), permuted = TRUE)
u_post = u_samp$u %>% data.frame() # (J * N_approx) x 2
u_mat = as.matrix(u_post)
u_post %>% dim
u_df = preprocess(u_post, D, prior)

hml_approx = hml_const(N_approx, D, u_df, J, prior)
hml_approx$const_vec


## bridge sampling stuff -------------------------------------------------------
log_density = function(u, data) {
    -psi(u, data)
}

lb <- rep(0, D)
ub <- rep(1, D)
colnames(u_mat) = names(u_df)[1:D]
names(lb) <- names(ub) <- colnames(u_mat)
## bridge sampling stuff -------------------------------------------------------
params = list(log_density = log_density, prior = prior, lb = lb, ub = ub)
bridge_result = bridgesampling::bridge_sampler(samples = u_mat,
                                               log_posterior = log_density,
                                               data = prior, lb = lb, ub = ub, silent = FALSE)
bridge_result$logml

hml_approx = hml_const(N_approx, D, u_df, J, prior)
hml_approx$const_vec # (N_approx x 1) vector


compute_bridge = function(samples, params) {
    bridge = bridgesampling::bridge_sampler(samples = samples,
                                            log_posterior = params$log_density,
                                            data = params$prior, 
                                            lb = params$lb, ub = params$ub, 
                                            silent = TRUE)
    bridge$logml
}




# STAN SETTINGS ----------------------------------------------------------------
J         = 2000         # number of MC samples per approximation
N_approx  = 100          # number of approximations
burn_in   = 2000         # number of burn in draws
n_chains  = 4            # number of markov chains to run
stan_seed = 123          # seed
J_iter = 1 / n_chains * N_approx * J + burn_in 
# STAN SETTINGS ----------------------------------------------------------------

# store approximations corresponding to each sample size
hyb = matrix(NA, N_approx, length(N_vec))
bridge = matrix(NA, N_approx, length(N_vec))

N_vec_log = seq(1, 15, by = 0.1)             # sample size grid unif in log
N_vec     = floor(exp(N_vec_log)) %>% unique # sample size to generate data
print(length(N_vec))               # number of different sample sizes


set.seed(1)
for (i in 105:length(N_vec)) {
    
    N = N_vec[i]   # pseudo-sample size
    
    gamma_dat = list(N = N) # for STAN sampler
    prior     = list(N = N) # for evaluation of psi, lambda
    
    # (1) generate posterior samples -- should give us (J * N_approx) draws
    gamma_fit_N = stan(file    =  stan_sampler, 
                       data    =  gamma_dat,
                       iter    =  J_iter,
                       warmup  =  burn_in,
                       chains  =  n_chains,                  
                       control =  list(adapt_delta = 0.99),  
                       refresh = 0)            
    
    u_samp = rstan::extract(gamma_fit_N, pars = c("u"), permuted = TRUE)
    u_post = u_samp$u %>% data.frame() # (J * N_approx) x 2
    
    #### (1) compute bridge estimate -------------------------------------------
    # prepare bridge sampling input
    u_mat = as.matrix(u_post)          # bridge_sampler() requires matrix input
    colnames(u_mat) = names(u_df)[1:D] # name columns of the samples
    
    # split u_mat into N_approx - many
    params = list(log_density = log_density, prior = prior, lb = lb, ub = ub)
    xsplit     = rep(1:N_approx, times = rep(J/N_approx, N_approx))
    u_mat_list = split.data.frame(u_mat, xsplit)
    bridge[,i] = unname(sapply(u_mat_list, compute_bridge, params = params))
    
    # mean(bridge[,i])
    #### (2) compute hybrid estimate -------------------------------------------
    u_df = preprocess(u_post, D, prior)
    hml_approx = hml_const(N_approx, D, u_df, J, prior)
    hyb[,i] = hml_approx$const_vec
    
    # (3) run algorithm to obtain N_approx estimates of the LIL
    print(paste("iter = ", i, "/", length(N_vec),  
                ' -- bridge = ', round(mean(bridge[,i]), 3),
                '; hyb = ', mean(hyb[,i]),
                sep = ''))
}



logZ = read.csv("singulard4.csv")


logn   = log(N_vec)
lil_df = data.frame(bridge = colMeans(bridge[,1:106]), 
                    hyb = colMeans(hyb[,1:106]), logn = logn[1:106],
                    logZ = logZ$logZ[1:106])

lil_df_long = melt(lil_df, id.vars = "logn")

formula1 = y ~ x

ggplot(lil_df_long, aes(x = logn, y = value, 
                        color = as.factor(variable))) + geom_point(size = 0.7) + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    labs(x = "log(n)", y = "log(Z)", 
         title = "Bridge (Red), Hybrid (Green), True (Blue)") + 
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16) + 
    theme(legend.position = "none")





