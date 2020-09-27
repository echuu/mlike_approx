

# setwd("C:/Users/ericc/Dropbox/logML") 
# source("setup.R")           # setup global environment, load in algo functions
# source("HIW_helper.R")

setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions
setwd("C:/Users/ericc/mlike_approx/covariance")
source("HIW_helper.R")      # covariance related helper functions
source("HIW_graph.R")



D = 9          # num rows/cols of the graph G and covariance/precision matrices
b = 3          # prior degrees of freedom
V = diag(1, D) # prior scale matrix 

D_0 = 0.5 * D * (D + 1) # num entries on diagonal and upper diagonal
J = 1e4


# Given graph testG
testG = matrix(c(1,1,0,0,1,0,0,0,0,
                 1,1,1,1,1,0,0,0,0,
                 0,1,1,1,0,0,0,0,0,
                 0,1,1,1,1,1,1,0,0,
                 1,1,0,1,1,1,0,0,0,
                 0,0,0,1,1,1,1,1,1,
                 0,0,0,1,0,1,1,1,1,
                 0,0,0,0,0,1,1,1,1,
                 0,0,0,0,0,1,1,1,1), D, D)

a = c(1, 3, 2, 5, 4, 6, 7, 8, 9)
testG = testG[a, a]


testG  = matrix(c(1,1,0,0,0,
                  1,1,1,1,1,
                  0,1,1,1,1,
                  0,1,1,1,1,
                  0,1,1,1,1), 5,5)


testG = matrix(c(1,1,1,0,0,
                  1,1,1,0,0,
                  1,1,1,1,1,
                  0,0,1,1,1,
                  0,0,1,1,1), 5,5)


D = nrow(testG)
b = 3          # prior degrees of freedom
V = diag(1, D) # prior scale matrix 

D_0 = 0.5 * D * (D + 1) # num entries on diagonal and upper diagonal
J = 1e4


# logical vector determining existence of edges between vertices
edgeInd = testG[upper.tri(testG, diag = TRUE)] %>% as.logical
upperInd = testG[upper.tri(testG)] %>% as.logical
D_u = sum(edgeInd)

# Specify true value of Sigma
set.seed(1)
true_params = HIWsim(testG, b, V)
Sigma_G = true_params$Sigma
Omega_G = true_params$Omega # precision matrix -- this is the one we work with

# chol(Omega_G)

# Generate data Y based on Sigma_G
N = 40
Y = matrix(0, N, D)
for (i in 1:N) {
    Y[i, ] = t(t(chol(Sigma_G)) %*% rnorm(D, 0, 1)) # (500 x D)
}

S = t(Y) %*% Y

params = list(N = N, D = D, D_0 = D_0, testG = testG, edgeInd = edgeInd,
              upperInd = upperInd, S = S, V = V, b = b)


## testing ---------------------------------------------------------------------

# calling sampleHIW() will return: post_samps, post_samps_0, Lt_post, 
# Sigma_post, Omega_post
postIW = sampleHIW(J, D_u, D_0, testG, b, N, V, S, edgeInd)  
post_samps = postIW$post_samps                 # (J x D_u)

# j = 37
# u = post_samps[j,] %>% unname %>% unlist
# u0 = postIW$post_samps_0[j,] %>% unname %>% unlist
# 
# HIW_logprior(u, params)  
# HIW_logprior_old(u0, params)
# 
# HIW_loglik(u, params)
# HIW_loglik_old(u0, params)
# 
# 
# 
# u0 = postIW$post_samps_0[1,] %>% unname %>% unlist
# # HIW_logprior_old(u0, params)
# 
# # reconstruct Lt_post
# postIW$Lt_post[[1]]
# 
# Lt_post_0 = matrix(0, D, D)
# Lt_vec_0 = numeric(D_0)
# Lt_vec_0[edgeInd] = u
# Lt_post_0[upper.tri(Lt_post_0, diag = T)] = Lt_vec_0
# Lt_post_0

# true log marginal likelihood
# logmarginal = function(Y, G, b, D, S){
#     n = nrow(Y); p = ncol(Y); # S = t(Y)%*%Y
#     logmarginal = -0.5*n*p*log(2*pi) + logHIWnorm(G, b, D) - 
#         logHIWnorm(G, b+n, D+S)
#     return(logmarginal)
# }


## testing ---------------------------------------------------------------------


u_df = preprocess(post_samps, D_u, params)     # J x (D_u + 1)

loglik_max = maxLogLik(Omega_G, params)

hml_approx = hml_const(1, D_u, u_df, J, params)

hml_approx$const_vec - loglik_max

# TODO: compute true log maximum likelihood
# lil(param_list) - loglik_max

# logmarginal(Y, testG, b, V) - loglik_max # using HIW exact

# TODO: compute log marginal likelihood with gnorm() function
# using gnorm in BDgraph to recalculate the same log Marginal likelihood
gnorm_approx = - 0.5 * D * N * log(2 * pi) + gnorm(testG, b + N, V + S, iter = 100) - 
    gnorm(testG, b, V, iter = 100)

# gnorm_wrong = -0.5 * D * N * log(2 * pi) + gnorm(testG, 3 + N, diag(1:D) + S, 100) - 
#     gnorm(testG, 3, diag(1:D), 100)


logmarginal(Y, testG, b, V, S)
gnorm_approx
hml_approx$const_vec




# hml_approx$param_out
# abs(logmarginal(Y, testG, b, V, S) - hml_approx$const_vec)
# abs(logmarginal(Y, testG, b, V, S) - gnorm_approx)
# hml_approx$param_out %>% 
#     dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs)

(orig_partition = hml_approx$param_out %>%
        dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs) %>% 
        dplyr::mutate(perc = n_obs / sum(n_obs)) %>% 
        arrange(desc(perc)) %>% 
        mutate(contrib = logQ_cstar / sum(logQ_cstar)))

K = nrow(orig_partition)

# ------------------------------------------------------------------------------

#### 2nd stage


## start here
n_samps = 10

# for the partition learned from prev fitted tree, extract the partition id and
# the optimal value of psi for this partition
og_part = hml_approx$param_out %>% 
    dplyr::select(-c(psi_choice, logQ_cstar))

# set.seed(1)
ss_part = fit_resid(og_part, D_u, n_samps, params)
ts_part = fit_resid(ss_part, D_u, n_samps / 2, params)
# fs_part = fit_resid(ts_part, D_u, n_samps / 2, params)


# truth

# original approx
x1 = hml_approx$const_vec
x2 = log_sum_exp(unlist(compute_expterms(ss_part, D_u)))
x3 = log_sum_exp(unlist(compute_expterms(ts_part, D_u)))
# log_sum_exp(unlist(compute_expterms(fs_part, D_u)))

mean(c(x1,x2,x3))
x2
x3
gnorm_approx
logmarginal(Y, testG, b, V, S)

log_sum_exp(unlist(compute_expterms(ts_part, D_u))) - logmarginal(Y, testG, b, V, S)


# ------------------------------------------------------------------------------
set.seed(2)
hml_approx = hml_const(1, D_u, u_df, J, params)
orig_partition = hml_approx$param_out %>%
        dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs) %>% 
        dplyr::mutate(perc = n_obs / sum(n_obs)) %>% 
        arrange(desc(perc)) %>% 
        mutate(contrib = logQ_cstar / sum(logQ_cstar))

K = nrow(orig_partition)
reapprox0 = resample_K(hml_approx, K, params, D_u)
hml_approx$const_vec
log_sum_exp(reapprox0$all_terms)
length(reapprox0$all_terms)

ts_approx_k = vector("list", K) 
ts_approx_terms = vector("list", K) 

for (k in 1:K) {
    
    print(paste("third stage on partition ", k, sep = ''))
    
    sub_part_k = reapprox0$ss_partitions[[k]]
    
    
    K_sub = nrow(sub_part_k$param_out)
    ts_approx_k[[k]] = resample_K(sub_part_k, K_sub, params, D_u, 5)
    
    ts_approx_terms[[k]] = ts_approx_k[[k]]$all_terms
}

log_sum_exp(unlist(ts_approx_terms))
length(unlist(ts_approx_terms))








# ------------------------------------------------------------------------------


J             = 2000                              # num MCMC samples from post
K_sims        = 50                                # num sims to run for each N
N_vec_log     = seq(5, 12, by = 0.25)             # sample size grid unif in log
N_vec         = floor(exp(N_vec_log)) %>% unique  # sample size to generate data
LIL_N_k_hat   = matrix(0, length(N_vec), K_sims)  # store hml() approximations
LIL_N_k_gnorm = matrix(0, length(N_vec), K_sims)  # store gnorm() approx
LIL_N         = numeric(length(N_vec))            # store true logML

length(N_vec)

lowerSigmaG = t(chol(Sigma_G))

i = 1; k = 1;

for (i in 1:length(N_vec)) {
    
    N = N_vec[i]
    
    LIL_N_k = numeric(K_sims) # store true log ML
    
    for (k in 1:K_sims) {
        
        Y = matrix(0, N, D)
        for (n in 1:N) {
            Y[n, ] = t(lowerSigmaG %*% rnorm(D, 0, 1)) # (500 x D)
        }
        
        S = t(Y) %*% Y
        
        params = list(N = N, D = D, D_0 = D_0, testG = testG, edgeInd = edgeInd,
                      upperInd = upperInd, S = S, V = V, b = b)
        
        # compute maximized log likelihood
        loglik_max = maxLogLik(Omega_G, params)
        
        # postIW contains: post_samps, Sigma_post, L_post
        postIW = sampleHIW(J, D_u, D_0, testG, b, N, V, S, edgeInd)  
        post_samps = postIW$post_samps                   # (J x D_u)
        
        u_df = preprocess(post_samps, D_u, params)       # J x (D_u + 1)
        
        # generate hybrid approximation
        # hml_approx() is a version of hml() that ignores the gradient term
        hml_approx = hml_const(1, D_u, u_df, J, params)
        
        # subtract maximized likelihood from the resulting approximation
        LIL_N_k_hat[i, k] = hml_approx$const_vec - loglik_max
        
        
        gnorm_approx = - 0.5 * D * N * log(2 * pi) + 
            gnorm(testG, b + N, V + S, iter = 100) - 
            gnorm(testG, b, V, iter = 100)
        
        LIL_N_k_gnorm[i, k] = gnorm_approx - loglik_max
        
        # hml_approx$const_vec - loglik_max # -34.70373
        
        # compute true log ML - maximized likelihood
        LIL_N_k[k] = logmarginal(Y, testG, b, V, S) - loglik_max
        
    }
    
    LIL_N[i] = mean(LIL_N_k)
    
    print(paste("iter ", i, "/", length(N_vec), ": ",
                "approx LIL for N = ", N, " -- LIL = ",
                round(mean(LIL_N_k_hat[i, ]), 2), 
                " (", round(LIL_N[i], 2), ")", 
                sep = ''))
    
}



library(reshape2)
library(ggpmisc)

formula1 = y ~ x

lil_0     = LIL_N                   # length(N_vec) x 1
lil_hyb   = rowMeans(LIL_N_k_hat)   # length(N_vec) x 1
lil_gnorm = rowMeans(LIL_N_k_gnorm) 
log_N     = log(N_vec)              # length(N_vec) x 1

LIL_df = data.frame(LIL_N = lil_0, LIL_hat = lil_hyb, LIL_g = lil_gnorm, 
                    log_N = log(N_vec))

LIL_df_long = melt(LIL_df, id.vars = "log_N")
head(LIL_df_long)

ggplot(LIL_df_long, aes(x = log_N, y = value, 
                        color = as.factor(variable))) + geom_point(size = 1.3) + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    labs(x = "log(n)", y = "log(Z)", 
         title = "True (Red), Approx (Blue), gnorm (color)") + 
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16) + 
    theme(legend.position = "none")












