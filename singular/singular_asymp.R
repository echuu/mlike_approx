


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

# path for dell
# DELL_PATH = "C:/Users/chuu/mlike_approx"
# setwd(DELL_PATH)

# source("partition/partition.R")         # load partition extraction functions
# source("hybrid_approx_v1.R")               # load main algorithm functions
# source("extractPartition.R")
source("C:/Users/ericc/mlike_approx/singular/singular_helper.R")    # load psi(), lambda()

# x11()

# STAN SETTINGS ----------------------------------------------------------------
J         = 5000         # number of MC samples per approximation
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

# (2) evaluate posterior samples using psi(u)
u_df_N = preprocess(u_post, D, prior)
hml_approx = hml_const(1, D, u_df_N, J, prior)
hml_approx$const_vec


plot(u_df_N[,1], u_df_N[,2], pch = 20, cex = 0.3, 
     col = rgb(0, 0, 0, alpha = 0.5),
     xlab = '', ylab = '', main = '')
rect(hml_approx$param_out$u1_lb, hml_approx$param_out$u2_lb,
     hml_approx$param_out$u1_ub, hml_approx$param_out$u2_ub, lwd = 3)



# x11()
ggplot(u_df_N, aes(u1, u2)) + geom_point()

# (3) run algorithm to obtain N_approx estimates of the LIL
# approx = hml(N_approx, D, u_df_N, J, prior)
# mean(approx) # -1.2044

# hml_approx = hml_const(1, D, u_df_N, J, prior)
hml_approx$const_vec

hml_approx$param_out %>%
    dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs)

# compute true value of logZ
n = N
library(pracma)
result = integral2(fun, 0, 1, 0, 1, reltol = 1e-50)
log(result$Q) # -1.223014 for n = 1000

# ------------------------------------------------------------------------------

# contains: leaf_id, u1_lb, u1_ub, ... , uD_lb, uD_ub, n_obs
part_0 = hml_approx$param_out %>% 
    dplyr::select(-c(psi_choice, psi_star, logQ_cstar))

part_set = part_0$leaf_id

(orig_partition = hml_approx$param_out %>%
        dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs) %>% 
        dplyr::mutate(perc = n_obs / sum(n_obs)))


K = length(part_set)

# initialize a list to store the vector containing the terms in exponential
# for each of the sub-partitions
# kth elmt is an s_k dim vector of terms that are to be exponentiated
# at the very end, all entries are unlisted and evaluated via log-sum-exp
exp_terms = vector("list", K) 
ck_star_list = vector("list", K)

perc_thresh = sort(orig_partition$perc, decreasing = T)

for (k in 1:K) {
    
    
    PERC_K = orig_partition[k,]$perc
    
    if (PERC_K >= perc_thresh[1]) {
        print("using original partition")
        # exp_terms[[k]] = hml_approx$const_approx[k]
        
        N_k_p = part_0$n_obs[k] * 10  # number of (re)samples to draw from part k
        part_k = part_0 %>%           # set of lower/upper bounds
            dplyr::filter(leaf_id == part_set[k]) %>%
            dplyr::select(-c(leaf_id, n_obs))
        
        # sample uniformly from each lower/upper bound pair to form a D-dim vector
        part_k_long = c(unlist(part_k)) %>% matrix(ncol = 2, byrow = T)
        
        resamp_k = Matrix_runif(N_k_p, lower = part_k_long[,1],
                                upper = part_k_long[,2]) %>% data.frame
        
        u_df_k = preprocess(resamp_k, D, prior) # N_k_p x (D_u + 1)
        
        c_k_approx = hml_const_mod(1, D, u_df_k, N_k_p, prior)
        
        ck_star_list[[k]] = c_k_approx$param_out %>%
            dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs)
        
        exp_terms[[k]] = c_k_approx$const_approx
        
    } else {
        N_k_p = part_0$n_obs[k] * 10  # number of (re)samples to draw from part k
        part_k = part_0 %>%           # set of lower/upper bounds
            dplyr::filter(leaf_id == part_set[k]) %>% 
            dplyr::select(-c(leaf_id, n_obs))
        
        # sample uniformly from each lower/upper bound pair to form a D-dim vector
        part_k_long = c(unlist(part_k)) %>% matrix(ncol = 2, byrow = T)
        
        resamp_k = Matrix_runif(N_k_p, lower = part_k_long[,1], 
                                upper = part_k_long[,2]) %>% data.frame
        
        u_df_k = preprocess(resamp_k, D, prior) # N_k_p x (D_u + 1)
        
        c_k_approx = hml_const(1, D, u_df_k, N_k_p, prior)
        
        ck_star_list[[k]] = c_k_approx$param_out %>%
            dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs) 
        
        exp_terms[[k]] = c_k_approx$const_approx
    }
    
}

all_terms = exp_terms %>% unlist

log_sum_exp(all_terms)    # -1.215615 (3/12), -1.230866 (2/12), -1.233479 (3/10)
hml_approx$const_vec      # -1.041315
log(result$Q)             # -1.223014 for n = 1000

abs(hml_approx$const_vec - lil(y, X, prior, post))
abs(log_sum_exp(all_terms) - log(result$Q))






























# library(tree)
# u_tree = tree(psi_u ~ ., u_df_N)
# plot(u_tree)
# text(u_tree, cex = 0.8)
# 
# plot(u_df_N[,1], u_df_N[,2], pch = 20, cex = 1, col = "cyan",
#      xlab = 'u1', ylab = 'u2', main = '')
# partition.tree(u_tree, add = TRUE, cex = 0.8, ordvars = c("u1", "u2"))


psi_fun = fun
singular_diag = approx_lil_diag(D, u_df_N, prior) 


singular_diag$logZ_numer
singular_diag$logZ_taylor1
singular_diag$lozZ_taylor2
singular_diag$verbose_partition

partition_info = singular_diag$partition_info %>% 
    mutate(numer = round(numer, 4), taylor1 = round(taylor1, 4), 
           lambda1 = round(lambda1, 5), lambda2 = round(lambda2, 5), 
          taylor2 = round(taylor2, 4), e_ck_2 = round(e_ck_2, 4))

write.csv(partition_info, "partition_info_singular.csv", 
          row.names = F)


plotPartition(u_df_N, singular_diag$verbose_partition)


# new stuff --------------------------------------------------------------------
singular_diag$logZ_numer
singular_diag$logZ_taylor1
singular_diag$lozZ_taylor2
singular_diag$hybrid

singular_diag$partition_approx
singular_diag$verbose_partition
singular_diag$param_out

singular_diag$n_const
singular_diag$n_taylor


source("hybrid_approx_v1.R")
test = hybrid_mlik(N_approx, D, u_df_N, J, prior)
test$const_vec  %>% mean
test$taylor_vec %>% mean
test$hybrid_vec %>% mean

# ------------------------------------------------------------------------------

# test generalized version of hybrid_mlik -- hml() function
hml_approx = hml(N_approx, D, u_df_N, J, prior)

# verify constant approximation
hml_approx$const_vec

# verify taylor approximation
hml_approx$taylor_vec

# verify hybrid approximation
hml_approx$hybrid_vec


rbind(hml_approx$const_vec,  hml_approx$const_vec_lse)
rbind(hml_approx$taylor_vec, hml_approx$taylor_vec_lse)
rbind(hml_approx$hybrid_vec, hml_approx$hybrid_vec_lse) # TODO




# run algorithm over grid of N -------------------------------------------------


# values of N for which we will compute + approximate the LIL
N_vec_log = seq(1, 17, by = 0.1)             # sample size grid unif in log
N_vec     = floor(exp(N_vec_log)) %>% unique # sample size to generate data

# N_vec     = c(24154953)         # sample size to use to generate data
print(length(N_vec))               # number of different sample sizes

# store approximations corresponding to each sample size
approx_taylor = matrix(NA, N_approx, length(N_vec))
approx_hybrid = matrix(NA, N_approx, length(N_vec))

set.seed(1)
for (i in 1:length(N_vec)) {
    
    N = N_vec[i]   # pseudo-sample size
    
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
    
    approx_out = hybrid_mlik(N_approx, D, u_df_N, J, prior) 
    
    approx_taylor[,i] = approx_out$taylor_vec
    approx_hybrid[,i] = approx_out$hybrid_vec
    
}


colMeans(approx_taylor)
colMeans(approx_hybrid)

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
logZ_taylor = colMeans(approx_taylor)  # 
logZ_hybrid = colMeans(approx_hybrid)  #
logn   = log(N_vec)

lil_df = data.frame(logZ_0 = logZ_0, logZ_taylor = logZ_taylor, 
                    logZ_hybrid = logZ_hybrid, logn = logn)

lil_df = lil_df[complete.cases(lil_df),]

lil_df_long = melt(lil_df, id.vars = "logn")




formula1 = y ~ x

ggplot(lil_df_long, aes(x = logn, y = value, 
                        color = as.factor(variable))) + geom_point(size = 0.7) + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    labs(x = "log(n)", y = "log(Z)", 
         title = "True (Red), Hybrid (Blue), Taylor (Green)") + 
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16) + 
    theme(legend.position = "none")

lm(logZ ~ logn, lil_df) # slope should be -0.25




# ------------------------------------------------------------------------------

# ggplot(u_df_N, aes(u1, u2)) + geom_point(col = 'red') + 
#     scale_x_continuous(name="x") + 
#     scale_y_continuous(name="y") +
#     annotate('rect', xmin = singular_diag$param_out$u1_lb[1:2], 
#              xmax = singular_diag$param_out$u1_ub[1:2], 
#              ymin = singular_diag$param_out$u2_lb[1:2], 
#              ymax = singular_diag$param_out$u2_ub[1:2], 
#              fill = alpha("grey", 0.2))
# 
# ggplot(partition_df) + 
#     geom_rect(aes(xmin = u1_lb, 
#                   xmax = u1_ub, 
#                   ymin = u2_lb, 
#                   ymax = u2_ub), alpha = 0.5,
#               fill = alpha("grey",0))







