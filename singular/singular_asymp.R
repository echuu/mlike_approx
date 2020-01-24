

library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio
library("ggpmisc")

# use stan to draw from the posterior distribution -----------------------------

# path for lenovo
LEN_PATH  = "C:/Users/ericc/mlike_approx"
setwd(LEN_PATH)

# path for dell
# DELL_PATH = "C:/Users/chuu/mlike_approx"
# setwd(DELL_PATH)

source("partition/partition.R")         # load partition extraction functions
source("hybrid_approx.R")               # load main algorithm functions
source("extractPartition.R")
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


# integral that we are trying to evaluate
fun <- function(x, y) exp(-n*x^2*y^4)


# one run of the algorithm -----------------------------------------------------
set.seed(1)
N = 500

gamma_dat = list(N = N) # for STAN sampler
prior     = list(N = N) # for evaluation of psi, lambda

# (1) generate posterior samples -- should give us (J * N_approx) draws
gamma_fit_N = stan(file    =  'singular/gamma_sample.stan', 
                   data    =  gamma_dat,
                   iter    =  J_iter,
                   warmup  =  burn_in,
                   chains  =  n_chains,                  
                   control =  list(adapt_delta = 0.99),  
                   refresh = 0)           

u_samp = rstan::extract(gamma_fit_N, pars = c("u"), permuted = TRUE)
u_post = u_samp$u %>% data.frame() # (J * N_approx) x 2

# (2) evaluate posterior samples using psi(u)
u_df_N = preprocess(u_post, D, prior)

ggplot(u_df_N, aes(u1, u2)) + geom_point()

# (3) run algorithm to obtain N_approx estimates of the LIL
approx = approx_lil(N_approx, D, u_df_N, J, prior)
mean(approx) # 

# compute true value of logZ
n = N
result = integral2(fun, 0, 1, 0, 1, reltol = 1e-50)
log(result$Q) # -1.223014 for n = 1000


# ------------------------------------------------------------------------------


# run algorithm over grid of N -------------------------------------------------


# values of N for which we will compute + approximate the LIL
N_vec_log = seq(1, 17, by = 0.1)             # sample size grid unif in log
N_vec     = floor(exp(N_vec_log)) %>% unique # sample size to generate data

# N_vec     = c(24154953)         # sample size to use to generate data
print(length(N_vec))               # number of different sample sizes

# store approximations corresponding to each sample size
approx_N = matrix(NA, N_approx, length(N_vec))

set.seed(1)
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
                       control =  list(adapt_delta = 0.99),  
                       refresh = 0)           
    
    u_samp = rstan::extract(gamma_fit_N, pars = c("u"), permuted = TRUE)
    u_post = u_samp$u %>% data.frame() # (J * N_approx) x 2
    
    # (2) evaluate posterior samples using psi(u)
    u_df_N = preprocess(u_post, D, prior)
    
    # (3) run algorithm to obtain N_approx estimates of the LIL
    
    print(paste("iter = ", i, "/", length(N_vec),  
                ' -- calculating logZ for N = ', N, 
                ' (', N_approx, ' approximations)', sep = ''))
    
    
    # N_approx, J settings indicate that J MC samples will use in each of the
    # N_approx estiamtes of the LIL -> return vector of length N_approx
    approx_N[,i] = approx_lil(N_approx, D, u_df_N, J, prior)
    
}

approx_N


colMeans(approx_N)

N_vec     = c(24154953)         # sample size to use to generate data
n = N_vec[1]
result = integral2(fun, 0, 1, 0, 1, reltol = 1e-50)
log(result$Q) # -1.223014 for n = 1000


# J = 1000 MC samples per approximation, N_approx = 10 approximations for 
# each sample size N -- grid below is needed to make use of the data
# N_vec_log = seq(1, 17, by = 0.1)             # sample size grid unif in log
# N_vec     = floor(exp(N_vec_log)) %>% unique # sample size to generate data
# write.csv(approx_N, "singular_asymptotics.csv", row.names = F)
# test_read = read.csv("singular_asymptotics.csv")

# ------------------------------------------------------------------------------

# numerical integration for the true log normalizing constant ------------------

library(pracma)

N_vec_log = seq(1, 17, by = 0.1)             # sample size grid unif in log
N_vec     = floor(exp(N_vec_log)) %>% unique # sample size to generate data

logZ_0 = rep(0,length(N_vec))

print(length(logZ_0))

i = 1;
for (n in N_vec) {
    result = integral2(fun, 0, 1, 0, 1, reltol = 1e-50)
    logZ_0[i] = log(result$Q)
    i = i + 1;
}

plot(logZ_0 ~ log(N_vec))
lm(logZ_0 ~ log(N_vec))


# ------------------------------------------------------------------------------

# overlay the true, approximate plotted vs. log n
library(reshape2)
logZ = colMeans(approx_N)  # approximation
logn   = log(N_vec)

lil_df = data.frame(logZ_0 = logZ_0, logZ = logZ, logn = logn)
lil_df_long = melt(lil_df, id.vars = "logn")


formula1 = y ~ x

ggplot(lil_df_long, aes(x = logn, y = value, 
                        color = as.factor(variable))) + geom_point(size = 0.7) + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    labs(x = "log(n)", y = "log(Z)", 
         title = "Approx (Blue), True Value via Numerical Integration (Red)") + 
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16) + 
    theme(legend.position = "none")

lm(logZ ~ logn, lil_df) # slope should be -0.25











