

# setup global environment, load in algo functions
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           

library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio
library("ggpmisc")
options(mc.cores = parallel::detectCores()) 

# use stan to draw from the posterior distribution -----------------------------

# path for lenovo
stan_sampler = 'C:/Users/ericc/mlike_approx/singular/gamma_sample.stan'
source("C:/Users/ericc/mlike_approx/singular/singular_helper.R")

# STAN SETTINGS ----------------------------------------------------------------
J         = 1e4          # number of MC samples per approximation
N_approx  = 1            # number of approximations
burn_in   = 2000         # number of burn in draws
n_chains  = 4            # number of markov chains to run
stan_seed = 123          # seed

J_iter = 1 / n_chains * N_approx * J + burn_in 
# ------------------------------------------------------------------------------

# GLOBAL MODEL SETTINGS --------------------------------------------------------
D = 2                    # dimension of parameter

# integral that we are trying to evaluate
fun <- function(x, y) exp(-n*x^2*y^4)

# one run of the algorithm -----------------------------------------------------
set.seed(123)
N = 1000

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
u_df = preprocess(u_post, D, prior)

hybrid = hybrid_ml(D, u_df, J, prior)
if (any(is.na(hybrid))) {print(paste("error in iteration", i)); next;}
hybrid$zhat

n = N
library(pracma)
result = integral2(fun, 0, 1, 0, 1, reltol = 1e-50)
log(result$Q)


N_vec_log = seq(6, 10, 0.02)        # sample size that is uniform over log scale
N_vec     = floor(exp(N_vec_log))   # sample size to use to generate data
logZ_0    = numeric(length(N_vec))  # store true value of log normalizing const
print(length(N_vec))                # number of different sample sizes


J = 1e4 * n_reps         # number of total MC samples
J_iter = 1 / n_chains * N_approx * J + burn_in 
n_reps = 100

approx_hybrid = matrix(NA, n_reps, length(N_vec))

for (i in 1:length(N_vec)) {
    
    N = N_vec[i]
    gamma_dat = list(N = N) # for STAN sampler
    prior     = list(N = N) # for evaluation of psi, lambda
    gamma_fit_N = stan(file    =  stan_sampler,
                       data    =  gamma_dat,
                       iter    =  J_iter,
                       warmup  =  burn_in,
                       chains  =  n_chains,
                       control =  list(adapt_delta = 0.99),
                       refresh = 0)
    u_samp = rstan::extract(gamma_fit_N, pars = c("u"), permuted = TRUE)
    u_post = u_samp$u %>% data.frame() # (J * N_approx) x 2
    u_df = preprocess(u_post, D, prior)
    # u_df %>% dim
    xsplit     = rep(1:n_reps, times = rep(J/n_reps, n_reps))
    u_df_list = split.data.frame(u_df, xsplit)
    
    for (j in 1:n_reps) {
        hybrid = hybrid_ml(D, u_df_list[[j]], J/n_reps, prior)
        if (any(is.na(hybrid))) {print(paste("error in iteration", i)); next;}
        approx_hybrid[j,i] = hybrid$zhat
    }
    
    n = N
    result = integral2(fun, 0, 1, 0, 1, reltol = 1e-50)
    logZ_0[i] = log(result$Q) 
    print(paste("iter = ", i, "/", length(N_vec), 
                " -- Calculating LIL for D = ", D, ", N = ", N, sep = ''))
    
} # end of loop iterating over different sample sizes

logZ_hybrid = colMeans(approx_hybrid)  #
logn = log(N_vec)


lil_df = data.frame(logZ_0 = logZ_0, logZ_hybrid = logZ_hybrid, logn = logn)

lil_df = data.frame(logZ_0 = logZ_0, logn = logn)

lil_df_long = melt(lil_df, id.vars = "logn")

saveRDS(list(J = J, lil_df = lil_df, approx_df = approx_hybrid), 
        file = 'singular_d2.RData')
singular_d2 = readRDS('singular_d2.RData')



formula1 = y ~ x

ggplot(lil_df_long, aes(x = logn, y = value, 
                        color = as.factor(variable))) + geom_point(size = 0.7) + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    labs(x = "log(n)", y = "log(Z)", 
         title = "True (Red), Hybrid (Blue), Taylor (Green), Constant (Purple)") + 
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16) + 
    theme(legend.position = "none")











