
library(rstan)
library(rstanarm)
library(ggplot2)
library(bayesplot)
library(dplyr)
library(gtools)    # for inv.logit() function

options(mc.cores = parallel::detectCores()) 

# save compiled stan program to hard disk so that don't need to recompile
rstan_options(auto_write = TRUE)

LEN_PATH  = "C:/Users/ericc/mlike_approx"
setwd(LEN_PATH)


source("partition/partition.R")
source("extractPartition.R")
source("hybrid_approx.R")
source("logistic/logistic_helper.R")

source("logistic/hme_logistic.R")


beta = sample(c(-5:5), D, replace = TRUE) # (0, 1, 3) for D = 3
# prior precision for beta
tau = 0.5
prior = list(y = y, X = X, D = D, N = N, tau = tau)

# number of MCMC samples
J = 4000

# ------------------------------------------------------------------------------

set.seed(1)

N = 300
D = 3
X = matrix(rnorm(N * D), nrow = N)
beta = sample(c(-5:5), D, replace = TRUE) # (0, 1, 3) for D = 3
z = X %*% beta
pr = exp(z) / (1 + exp(z))
y = rbinom(N, 1, pr)

df = data.frame(y, X)

# prior precision for beta
tau = 0.5
prior = list(y = y, X = X, D = D, N = N, tau = tau)

# number of MCMC samples
J = 4000

## TODO: figure out how to specify the number of posterior samples used
bglm_fit = stan_glm(y ~ . -1, data = df, 
                    prior = normal(location = 0, scale = 1 / sqrt(tau)),
                    family = binomial(link = 'logit'),
                    refresh = 0)
summary(bglm_fit)

u_samps = bglm_fit %>% as.data.frame() # extract posterior samples

psi_u = apply(u_samps, 1, psi, prior = prior)
u_df = preprocess(u_samps, D, prior)

hml_approx = hml(1, D, u_df, J, prior) 

loglik_J = logistic_lik(u_df, prior, J, D, N)


hist(hml_approx$taylor_vec - loglik_J)

hml_approx$const_vec
hml_approx$taylor_vec
hml_approx$hybrid_vec
hml_approx$partition


# ------------------------------------------------------------------------------
set.seed(1)
K_sims = 10            # num of simulations to run FOR EACH N in N_vec
N_vec_log = seq(5, 10, by = 0.25)            # sample size grid unif in log
N_vec     = floor(exp(N_vec_log)) %>% unique # sample size to generate data

rbind(N_vec_log, N_vec)


# K_sims = 1           # num of simulations to run FOR EACH N in N_vec
# N_vec = c(100)
LIL_N_k = matrix(0, length(N_vec), K_sims) # store the K_sims results

for (i in 1:length(N_vec)) {
    
    N = N_vec[i]
    
    for (k in 1:K_sims) {
        
        # generate data --------------------------------------------------------
        X = matrix(rnorm(N * D), nrow = N)
        z = X %*% beta
        pr = exp(z) / (1 + exp(z))
        y = rbinom(N, 1, pr)
        
        
        prior = list(y = y, X = X, D = D, N = N, tau = tau)
        
        df = data.frame(y, X)
        
        # obtain posterior samples ---------------------------------------------
        bglm_fit = stan_glm(y ~ . -1, data = df, 
                            prior = normal(location = 0, scale = 1 / sqrt(tau)),
                            family = binomial(link = 'logit'),
                            refresh = 0)
        
        u_samps = bglm_fit %>% as.data.frame() # extract posterior samples
        u_df = preprocess(u_samps, D, prior)
        
        # run algorithm --------------------------------------------------------
        
        hml_k = hml(1, D, u_df, J, prior)
        
        # subtract the loglikelihood
        LIL_N_k[i, k] = hml_k$taylor_vec -
            (sum(y * (X %*% beta) - log(1 + exp(X %*% beta))))
    }
    
    print(paste("approximating LIL for N = ", N, " -- LIL = ", 
                round(mean(LIL_N_k[i, k]), 2), sep = ''))
    
    
}


apply(LIL_N_k, 1, mean, na.rm = T)
apply(LIL_N_k, 1, var, na.rm = T)

LIL_df = data.frame(LIL_N = rowMeans(LIL_N_k), log_N = log(N_vec))
plot(LIL_N ~ log_N, LIL_df)

lm(LIL_N ~ log_N, LIL_df)


library(reshape2)
library(ggpmisc)
formula1 = y ~ x

ggplot(LIL_df, aes(x = log_N, y = LIL_N)) + geom_point(size = 1.5) + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    labs(x = "log(n)", y = "log(ML)", 
         title = "Approximate ML for Logistic Regression, D = 3") + 
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16) + 
    theme(legend.position = "none")











# LIL_const[,k]  = hml_approx$const_vec
# LIL_taylor[,k] = hml_approx$taylor_vec
# LIL_hybrid[,k] = hml_approx$hybrid_vec



# https://stats.stackexchange.com/questions/209810/computation-of-the-marginal-likelihood-from-mcmc-samples

# https://arxiv.org/pdf/0907.5123.pdf

# https://warwick.ac.uk/fac/sci/statistics/crism/workshops/estimatingconstants

# https://www.cse.iitk.ac.in/users/piyush/courses/tpmi_winter19/tpmi_w19_lec7_slides_print.pdf

