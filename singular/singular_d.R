
# setup global environment, load in algo functions
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           
source("C:/Users/ericc/mlike_approx/singular/singular_helper.R")

library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio
library("ggpmisc")
options(mc.cores = parallel::detectCores()) 

# path for lenovo
stan_sampler = 'C:/Users/ericc/mlike_approx/singular/gamma_sample.stan'

# STAN SETTINGS ----------------------------------------------------------------
J         = 2000         # number of MC samples per approximation
N_approx  = 1            # number of approximations
burn_in   = 2000         # number of burn in draws
n_chains  = 4            # number of markov chains to run
stan_seed = 123          # seed

J_iter = 1 / n_chains * N_approx * J + burn_in 

D = 2

# one run of the algorithm -----------------------------------------------------
set.seed(123)
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

## bridge sampling stuff -------------------------------------------------------
log_density = function(u, data) {
    -psi(u, data)
}

lb <- rep(0, D)
ub <- rep(1, D)
colnames(u_mat) = names(u_df)[1:D]
names(lb) <- names(ub) <- colnames(u_mat)
params = list(log_density = log_density, prior = prior, lb = lb, ub = ub)
## bridge sampling stuff -------------------------------------------------------

bridge_result = bridgesampling::bridge_sampler(samples = u_mat,
                               log_posterior = log_density,
                               data = prior, lb = lb, ub = ub, silent = FALSE)
bridge_result$logml

compute_bridge = function(samples, params) {
    bridge = bridgesampling::bridge_sampler(samples = samples,
                                            log_posterior = params$log_density,
                                            data = params$prior, 
                                            lb = params$lb, ub = params$ub, 
                                            silent = TRUE)
    bridge$logml
}


## numerical integration -------------------------------------------------------
n = N
library(pracma)
fun <- function(x, y) exp(-n*x^2*y^4)
result = integral2(fun, 0, 1, 0, 1, reltol = 1e-50)
log(result$Q) # -1.223014 for n = 1000

hml_approx = hml_const(N_approx, D, u_df, J, prior)
hml_approx$const_vec # (N_approx x 1) vector
# print(paste(nrow(u_df) / N_approx, " MCMC samples per approximation", sep = ''))

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

N_vec_log = seq(1, 9, by = 0.1)             # sample size grid unif in log
N_vec     = floor(exp(N_vec_log)) %>% unique # sample size to generate data
print(length(N_vec))               # number of different sample sizes


set.seed(1)
for (i in 1:length(N_vec)) {
    
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
    # u_df = preprocess(u_post, D, prior)
    # hml_approx = hml_const(N_approx, D, u_df, J, prior) 
    # hyb[,i] = hml_approx$const_vec
    
    # (3) run algorithm to obtain N_approx estimates of the LIL
    print(paste("iter = ", i, "/", length(N_vec),  
                ' -- bridge = ', round(mean(bridge[,i]), 3),
                # '; hyb = ', round(mean(hyb[,i]), 3),
                ' (', N_approx, ' approximations)', sep = ''))
}



#### (4) plot results
library(pracma)
N_vec_log = seq(1, 9, by = 0.1)             # sample size grid unif in log
N_vec     = floor(exp(N_vec_log)) %>% unique # sample size to generate data
logZ_0 = rep(0,length(N_vec))
print(length(logZ_0))
fun <- function(x, y) exp(-n*x^2*y^4)
i = 1;
for (n in N_vec) {
    result = integral2(fun, 0, 1, 0, 1, reltol = 1e-50)
    logZ_0[i] = log(result$Q)
    i = i + 1;
}

plot(logZ_0 ~ log(N_vec))
lm(logZ_0 ~ log(N_vec))
logn   = log(N_vec)
lil_df = data.frame(logZ_0 = logZ_0, bridge = colMeans(bridge), logn = logn)

lil_df_long = melt(lil_df, id.vars = "logn")

formula1 = y ~ x

ggplot(lil_df_long, aes(x = logn, y = value, 
                        color = as.factor(variable))) + geom_point(size = 0.7) + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    labs(x = "log(n)", y = "log(Z)", 
         title = "True (Red), bridge (Blue)") + 
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16) + 
    theme(legend.position = "none")


df_out = data.frame(n = N_vec, logn = logn, logZ_0 = logZ_0, 
                    bridge = colMeans(bridge))

write.csv(df_out, "bridge_d2.csv", row.names = F)
colnames(bridge) = paste("n_", N_vec, sep = '')
write.csv(bridge, "bridge_d2_raw.csv", row.names = F)
read.csv("bridge_d2_raw.csv") %>% head


# 
# xmat <- matrix(1:100, nrow=20, ncol=5, byrow=TRUE)
# xsplit <- rep( 1:5, times= rep(4,5))
# tmp <- split.data.frame(xmat,xsplit)
# 
# N_approx = 10
# xsplit <- rep(1:N_approx, times = rep(J/N_approx, N_approx))
# x_mat_list <- split.data.frame(u_mat, xsplit)
# 
# xsplit     = rep(1:N_approx, times = rep(J/N_approx, N_approx))
# u_df_list = split.data.frame(u_df, xsplit)
#     
# set.seed(1)
# unlist(unname(sapply(x_mat_list, compute_bridge, params = params)))
# 
# unname(sapply(u_df_list, hml_const_vec, N_approx = N_approx, D = D, 
#               J = J/N_approx, prior = prior))


