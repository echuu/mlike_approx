
# covarIW.R 


# change to your working directory
# note: all files must be in the same directory
# setwd("C:/Users/ericc/Dropbox/logML")

# important: files below must be sourced in order
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions
setwd("C:/Users/ericc/mlike_approx/covariance")
source("covarIW_helper.R")  # covariance related helper functions


N = 100                     # number of observations
D = 10                      # num rows/cols in the covariance matrix
D_u = 0.5 * D * (D + 1)     # dimension of u that is fed into the tree
J = 2000


## wishart prior parameters
Omega = diag(1, D)          # scale matrix
nu    = D + 1               # degrees of freedom

## specify the true covariance matrix
set.seed(1)
Sigma = matrix(rWishart(1, D, Omega), D)
is.positive.definite(Sigma)


##
## below is one iteration of the simulation: 
##     (1) generate data
##     (2) sample from posterior
##     (3) compute: 
##                  (a) maximized likelihood 
##                  (b) approximate logML
##                  (c) true logML
## 

## (1) generate data
X = rmvnorm(N, mean = rep(0, D), sigma = Sigma) # (N x p)
S = t(X) %*% X                                  # (p x p)


## store parameters in a list that can be passed into the algorithm
param_list = list(S = S, N = N, D = D, D_u = D_u, # S, dimension vars
                  Omega = Omega, nu = nu)         # prior params


## (2) obtain posterior samples
postIW = sampleIW(J, N, D_u, nu, S, Omega)     # post_samps, Sigma_post, L_post


# these are the posterior samples stored as vectors (the lower cholesky factors
# have been collapsed into D_u dimensional vectors)
post_samps = postIW$post_samps                 # (J x D_u)


# u_df stores the posterior samples row-wise so that the first D_u columns 
# store the lower cholesky factors in vector form, and the last column is
# the function evaluate psi(u), so u \in R^(D_u), and psi(u) \in R
u_df = preprocess(post_samps, D_u, param_list) # J x (D_u + 1)


## (3a) compute maximized likelihood (using true Sigma)
loglik_max = maxLogLik(Sigma, param_list)


## (3b) compute approximation
hml_approx = hml_const(1, D_u, u_df, J, param_list)

hml_approx$param_out %>%
    dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs)


# the log ML approximation is stored in the "const_vec" variable
# subtract off the maximized log likelihood
# hml_approx$const_vec - loglik_max


# (3c) compute true log ML, subtract off maximized log likelihood
# lil(param_list) - maxLogLik(Sigma, param_list)

(approx_logml = hml_approx$const_vec) # -1178.765
(true_logml = lil(param_list))        # -1187.283

# abs(approx_logml - true_logml)

# 2nd stage sampling -----------------------------------------------------------


# (orig_partition = hml_approx$param_out %>%
#      dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs) %>% 
#      dplyr::mutate(perc = n_obs / sum(n_obs)) %>% 
#      arrange(desc(perc)) %>% 
#      mutate(contrib = logQ_cstar / sum(logQ_cstar)))

K = nrow(orig_partition)

# set.seed(1)

hml_obj = hml_approx
hml_obj = NULL
n_samps = 10

reapprox0 = resample_K(hml_approx, K, param_list, D_u)
log_sum_exp(reapprox0$all_terms)
# length(reapprox0$all_terms)


# 3rd stage sampling -----------------------------------------------------------











# 2nd/3rd stage sampling w/ reps -----------------------------------------------



orig_approx = numeric(G)
ss_approx   = numeric(G)
ts_approx   = numeric(G)
true_logml  = numeric(G)
for (g in 1:G) {
    
    X = rmvnorm(N, mean = rep(0, D), sigma = Sigma) # (N x p)
    S = t(X) %*% X                                  # (p x p)
    
    
    ## store parameters in a list that can be passed into the algorithm
    param_list = list(S = S, N = N, D = D, D_u = D_u, # S, dimension vars
                      Omega = Omega, nu = nu)         # prior params
    
    
    true_logml[g] = lil(param_list)
    
    ## (2) obtain posterior samples
    postIW = sampleIW(J, N, D_u, nu, S, Omega)     # post_samps, Sigma_post, L_post
    
    
    # these are the posterior samples stored as vectors (the lower cholesky factors
    # have been collapsed into D_u dimensional vectors)
    post_samps = postIW$post_samps                 # (J x D_u)
    
    
    # u_df stores the posterior samples row-wise so that the first D_u columns 
    # store the lower cholesky factors in vector form, and the last column is
    # the function evaluate psi(u), so u \in R^(D_u), and psi(u) \in R
    u_df = preprocess(post_samps, D_u, param_list) # J x (D_u + 1)
    
    
    ## (3a) compute maximized likelihood (using true Sigma)
    loglik_max = maxLogLik(Sigma, param_list)
    
    
    ## (3b) compute approximation
    hml_approx = hml_const(1, D_u, u_df, J, param_list)
    
    orig_partition = hml_approx$param_out %>%
        dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs) %>%
        dplyr::mutate(perc = n_obs / sum(n_obs)) %>%
        arrange(desc(perc)) %>%
        mutate(contrib = logQ_cstar / sum(logQ_cstar))
    
    K = nrow(orig_partition)
    
    orig_approx[g] = hml_approx$const_vec 
    
    reapprox0 = resample_K(hml_approx, K, param_list, D_u)
    ss_approx[g] = log_sum_exp(reapprox0$all_terms)
    
    
    ts_approx_k = vector("list", K) 
    ts_approx_terms = vector("list", K) 
    
    for (k in 1:K) {
        
        # print(paste("third stage on partition ", k, sep = ''))
        
        sub_part_k = reapprox0$ss_partitions[[k]]
        
        K_sub = nrow(sub_part_k$param_out)
        ts_approx_k[[k]] = resample_K(sub_part_k, K_sub, param_list, D_u, 5)
        
        ts_approx_terms[[k]] = ts_approx_k[[k]]$all_terms
    }
    
    ts_approx[g] = log_sum_exp(unlist(ts_approx_terms))
    
    
    print(paste('iter ', g, '/', G, ' : ', 
                round(ss_approx[g], 4), ' (err: ', 
                round(abs(true_logml[g] - ss_approx[g]), 4), ', avg: ', 
                round(mean(ss_approx[1:g]), 4), '), ',
                round(ts_approx[g], 4), ' (err: ', 
                round(abs(true_logml[g] - ts_approx[g]), 4), ', avg: ', 
                round(mean(ts_approx[1:g]), 4), '), ', sep = ''))
    
}

mean(orig_approx) # no re-partioning
mean(ss_approx)   # second stage 
mean(ts_approx)   # third stage
mean(true_logml)  # true logML


abs(mean(true_logml) - mean(orig_approx))
abs(mean(true_logml) - mean(ss_approx))
abs(mean(true_logml) - mean(ts_approx))




### asymptotics ----------------------------------------------------------------

# ------------------------------------------------------------------------------
J           = 1e5                               # num of MCMC samples from post
K_sims      = 25                                # num of sims to run for each N
N_vec_log   = seq(5, 12, by = 0.25)             # sample size grid unif in log
N_vec       = floor(exp(N_vec_log)) %>% unique  # sample size to generate data
LIL_N_k_hat = matrix(0, length(N_vec), K_sims)  # store approximations
LIL_N       = numeric(length(N_vec))            # store true logML

length(N_vec)

for (i in 1:length(N_vec)) {
    
    N = N_vec[i]
    
    LIL_N_k = numeric(K_sims) # store true log ML
    
    for (k in 1:K_sims) {
        # set.seed(1)
        X = rmvnorm(N, mean = rep(0, D), sigma = Sigma) # (N x p)
        S = matrix(0, D, D)
        for (n in 1:N) {
            S = S + tcrossprod(X[n,]) # compute sum_n x_n * x_n'
        }
        
        param_list = list(S = S, N = N, D = D, D_u = D_u, # S, dimension vars
                          Omega = Omega, nu = nu,         # prior parameters
                          u_df = NULL)                    # posterior samples
        
        # compute maximized log likelihood
        loglik_max = maxLogLik(Sigma, param_list)
        
        # postIW contains: post_samps, Sigma_post, L_post
        postIW = sampleIW(J, N, D_u, nu, S, Omega) 
        
        post_samps = postIW$post_samps                   # (J x D_u)
        u_df = preprocess(post_samps, D_u, param_list)   # J x (D_u + 1)
        
        # generate hybrid approximation
        # hml_approx() is a version of hml() that ignores the gradient term
        hml_approx = hml_const(1, D_u, u_df, J, param_list)
        
        # subtract maximized likelihood from the resulting approximation
        LIL_N_k_hat[i, k] = hml_approx$const_vec - loglik_max
        
        # hml_approx$const_vec - loglik_max # -34.70373
        
        # compute true log ML - maximized likelihood
        LIL_N_k[k] = lil(param_list) - loglik_max
        
        # lil(param_list) - loglik_max # -43.2191
        
    }
    
    LIL_N[i] = mean(LIL_N_k)
    
    print(paste("iter ", i, "/", length(N_vec), ": ",
                "approx LIL for N = ", N, " -- LIL = ",
                round(mean(LIL_N_k_hat[i, ]), 2), 
                " (", round(LIL_N[i], 2), ")", 
                sep = ''))
    
}



lil_hyb   = rowMeans(LIL_N_k_hat)  # length(N_vec) x 1
log_N     = log(N_vec)             # length(N_vec) x 1

LIL_df = data.frame(LIL_hat = lil_hyb, log_N = log_N)
# LIL_df = data.frame(LIL_N = LIL_N, log_N = log(N_vec))
write.csv(LIL_df, "covarIW_J1e5.csv", row.names = F) # N = seq(5, 12, by = 0.25)



library(reshape2)
library(ggpmisc)

formula1 = y ~ x


# approx
ggplot(LIL_df, aes(x = log_N, y = LIL_hat)) + geom_point(size = 1.3) + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    labs(x = "log(n)", y = "log(Z)", 
         title = "Hybrid (Blue)") + 
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8)


# true
ggplot(LIL_df, aes(x = log_N, y = LIL_N)) + geom_point(size = 1.3) + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    labs(x = "log(n)", y = "log(Z)", 
         title = "Hybrid (Blue)") + 
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8)


# TODO: with the true log marginal likelihood
lil_0   = LIL_N                  # true value,   length(N_vec) x 1
lil_hyb = rowMeans(LIL_N_k_hat)  # approx value, length(N_vec) x 1
log_N   = log(N_vec)             # log(N) grid,  length(N_vec) x 1

LIL_df = data.frame(LIL_N = lil_0, LIL_hat = lil_hyb, log_N = log(N_vec))

LIL_df_long = melt(LIL_df, id.vars = "log_N")
head(LIL_df_long)

ggplot(LIL_df_long, aes(x = log_N, y = value, 
                        color = as.factor(variable))) + geom_point(size = 1.3) + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    labs(x = "log(n)", y = "log(Z)", 
         title = "True (Red), Approx (Blue)") + 
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16) + 
    theme(legend.position = "none")
































#### model output diagnostics

part_info = hml_approx$param_out
u_df_info = hml_approx$u_df_fit

# (1.1) for each point: fitted value, leaf_id, residual

## notation: 
# psi_u     : true value, psi(u)
# psi_star  : psi(u_star), where u_star is the partition's representative point
# psi_resid : psi_u - psi_star
# psi_hat   : fitted value for psi(u) on a given partition
##

u_df_info %>% head

# (1.2) for each partition: 'median' points, fitted value, upper/lower bounds,
# number of observations in partition
part_info 


# compare the psi values from tree vs. psi values evaluated at 'representative'
# point of each partition
part_psi = part_info %>% 
    dplyr::select(leaf_id, psi_hat) %>% 
    merge(u_df_info %>% dplyr::select(leaf_id, psi_star) %>% unique, 
          by = 'leaf_id')

part_psi

# (2) fitted value for each partition (const_approx)
# (2.1) look at the MSE for each of the partitions to see if there's one
# parition where there is significantly larger discrepancy
# (2.2) consider the number of observations in each of the partitions in
# conjunction with the MSE for each partition

library(MLmetrics)
part_mse = u_df_info %>% 
    dplyr::group_by(leaf_id) %>% 
    summarise(psi_star_mse = MSE(psi_u, psi_star)) %>% 
    merge(part_info %>% dplyr::select(leaf_id, n_obs), by = 'leaf_id')


# for each partition, display the fitted value for psi (from tree), 
# psi evaluated at the representative point, the mse associated with the
# representative point, the number of observations for that partition
psi_df = merge(part_psi, part_mse, by = 'leaf_id')

psi_df %>% arrange(psi_star_mse)

u_df_info %>% filter(leaf_id == 28)

psi_df %>% arrange(n_obs)

# compute mse for the tree's fitted values
psi_rpart = hml_approx$u_rpart
tree_dev = rpart_mse(psi_rpart, u_df)
tree_mse = tree_dev %>% group_by(leaf_id) %>% 
    summarise(psi_hat_mse = mean(dev_sq))

# gather all psi approximations (tree and algo) with associated MSEs
# into one table
psi_mse = psi_df %>% merge(tree_mse, by = 'leaf_id') %>% 
    dplyr::select(leaf_id, psi_hat, psi_hat_mse, psi_star, psi_star_mse, n_obs)


psi_mse %>% arrange(n_obs) %>% mutate(n_perc = n_obs / J)


hml_approx$const_vec
lil(param_list)



# (3) compute approximate integral over each partition
hml_approx$const_approx


# (4) extract partitions corresponding to the diagonal

# extract column names corresponding to diagonal entries (this is built into 
# getDiagCols() function already, but could be of some use later to have these
# indices available
diagInd = getDiagIndex(D, D_u)
diag_names = paste("u", diagInd, sep = '')

# extract u_k_star with corresponding lb, ub for each k that is a diagonal entry
getDiagCols(part_info, D, D_u)













# ------------------------------------------------------------------------------

library(cubature)
library(matrixcalc)


e_psi = function(u) {
    
    exp(-psi(u, param_list))
    
}

u_0 = u_df[1,-4] %>% unname %>% unlist
L_test = matrix(0, D, D)
L_test[lower.tri(L_test, diag = T)] = u_0

stable_numer = adaptIntegrate(e_psi, 
                              lowerLimit = c(0, -Inf, -Inf, 0, -Inf, 0), 
                              upperLimit = rep(Inf, D_u),
                              tol = 1e-4)

hml_approx$const_vec
log(stable_numer$integral)
lil(param_list)











