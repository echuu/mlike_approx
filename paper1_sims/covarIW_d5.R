
# important: files below must be sourced in order
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions
setwd("C:/Users/ericc/mlike_approx/covariance")
source("covarIW_helper.R")  # covariance related helper functions


source("C:/Users/ericc/mlike_approx/covariance/covarIW_helper.R")
Rcpp::sourceCpp("C:/Users/ericc/mlike_approx/speedup/fast_covIW.cpp")

N = 100                     # number of observations
D = 6                       # num rows/cols in the covariance matrix
D_u = 0.5 * D * (D + 1)     # dimension of u that is fed into the tree
J = 1000


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
X = mvtnorm::rmvnorm(N, mean = rep(0, D), sigma = Sigma) # (N x p)
S = t(X) %*% X                                  # (p x p)


## store parameters in a list that can be passed into the algorithm
param_list = list(S = S, N = N, D = D, D_u = D_u, # S, dimension vars
                  Omega = Omega, nu = nu)         # prior params

# ------------------------------------------------------------------------------

log_density = function(u, data) {
    -psi(u, data)
}

hme_exp_term = function(u) {

    L = matrix(0, D, D)             # (D x D) lower triangular matrix
    L[lower.tri(L, diag = T)] = u   # populate lower triangular terms

    logDiagL = log(diag(L))         # log of diagonal terms of L


    # compute loglikelihood using the 'L' instead of 'Sigma' --> allows us
    # to use simplified version of the determinant

    loglik = - 0.5 * N * D * log(2 * pi) - N * sum(logDiagL) -
        0.5 * matrix.trace(solve(L %*% t(L)) %*% S)

    -loglik
}


hme_approx = function(u_df) {

    # N   = params$N
    # D   = params$D
    # D_u = params$D_u

    tmp_J = apply(u_df[,1:D_u], 1, hme_exp_term) # term in the exponential
    log(J) - log_sum_exp(tmp_J)
}





# ------------------------------------------------------------------------------

## (2) obtain posterior samples
postIW = sampleIW(J, N, D_u, nu, S, Omega)     # post_samps, Sigma_post, L_post

# these are the posterior samples stored as vectors (the lower cholesky factors
# have been collapsed into D_u dimensional vectors)
post_samps = postIW$post_samps                 # (J x D_u)
u_df = preprocess(post_samps, D_u, param_list)


u_samp = as.matrix(post_samps)
colnames(u_samp) = names(u_df)[1:D_u]
# prepare bridge_sampler input()
lb = rep(-Inf, D_u)
ub = rep(Inf, D_u)

# diag_ind = getDiagIndex(D, D_u) # obtain column index of the diagonal entries
# lb[diag_ind] = 0                # diagonal entries are positive
names(lb) <- names(ub) <- colnames(u_samp)

bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
                                               log_posterior = log_density,
                                               data = param_list,
                                               lb = lb, ub = ub,
                                               silent = TRUE)
bridge_result$logml



(LIL = lil(param_list))
# postIW = sampleIW(J, N, D_u, nu, S, Omega)     # post_samps, Sigma_post, L_post
# post_samps = postIW$post_samps                 # (J x D_u)
# u_df = preprocess(post_samps, D_u, param_list) # J x (D_u + 1)
hybrid = hybrid_ml(D_u, u_df, J, param_list)
hybrid$zhat

lambda = function(u, params) { pracma::grad(psi_covar, u, params = params) }
hess   = function(u, params) { pracma::hessian(psi_covar, u, params = params) }

hybridml::hybml_const(u_df)$zhat
hybridml::hybml(u_df, param_list, grad = lambda, hess = hess)


# try the updated approximation using EP
hyb_numer(u_df, psi = psi_covar, params = param_list)
hyb(u_df, psi = psi_covar, params = param_list)




hme_approx(u_df)


B = 100 # number of replications
hyb    = numeric(B)
bridge = numeric(B)
hme    = numeric(B)
J = 25
set.seed(1)
i = 1
while (i <= B) {

    postIW = sampleIW(J, N, D_u, nu, S, Omega)     # post_samps, Sigma_post, L_post

    # these are the posterior samples stored as vectors (the lower cholesky factors
    # have been collapsed into D_u dimensional vectors)
    post_samps = postIW$post_samps                 # (J x D_u)


    # u_df stores the posterior samples row-wise so that the first D_u columns
    # store the lower cholesky factors in vector form, and the last column is
    # the function evaluate psi(u), so u \in R^(D_u), and psi(u) \in R
    u_df = preprocess(post_samps, D_u, param_list) # J x (D_u + 1)
    # u_df_fast = preprocess(post_samps, D_u, param_list) # J x (D_u + 1)
    #
    # ## (3b) compute approximation
    # hybrid = logml(D_u, u_df, J, param_list)
    # hybrid = hybrid_ml(D_u, u_df, J, param_list)
    # if (any(is.na(hybrid))) {print(paste("error in iteration", i)); next;}
    #
    # hyb[i] = hybrid$zhat


    u_samp = as.matrix(post_samps)
    colnames(u_samp) = names(u_df)[1:D_u]
    lb = rep(-Inf, D_u)
    ub = rep(Inf, D_u)
    names(lb) <- names(ub) <- colnames(u_samp)

    bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
                                                   log_posterior = log_density,
                                                   data = param_list,
                                                   lb = lb, ub = ub,
                                                   silent = TRUE,
                                                   method = 'warp3')
    bridge[i] = bridge_result$logml
    # hme[i] = hme_approx(u_df)

    print(paste("iter ", i, ': ',
                # "hybrid = ", round(mean(hyb[hyb!=0]), 3),
                # '; ', "ae = ", round(mean((LIL - hyb[hyb!=0])), 4),
                # ' // ',
                "bridge = ", round(mean(bridge[bridge!=0]), 3), '; ',
                "ae = ", round(mean((LIL - bridge[bridge!=0])), 4),
                # "hme = ", round(mean(hme[hme!=0]), 3), '; ',
                # "ae = ", round(mean((LIL - hme[hme!=0])), 4),
                sep = ''))
    i = i + 1

}


approx = data.frame(LIL, hyb = hyb, bridge = bridge, hme = hme)
error = data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
                   ae = colMeans((LIL - approx)),
                   rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)
error

fig_loc = 'C:/Users/ericc/mlike_approx/final_sims/'
saveRDS(list(J = J, D = D, D_u = D_u, N = N, approx = approx, error = error),
        file = paste(fig_loc, 'covarIW_d4_j25.RData', sep = ''))

covar_iw_d4 = readRDS(paste(fig_loc, 'covarIW_d4_j25.RData', sep = ''))

approx = covar_iw_d4$approx
approx %>% head
approx$wbse = bridge

names(approx) = c("LIL", "HybE", "BSE", "HME", "WBSE")

approx = approx %>% dplyr::select("LIL", "HybE", "BSE", "WBSE", "HME")

(error = data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
                    ae = colMeans(LIL - approx),
                    rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3))

delta_df = (LIL - approx) %>% dplyr::select(-c('LIL')) %>% melt()

delta_df %>% head

### figure dimensions for saving: (954 x 488)
ggplot(delta_df, aes(x = variable, y = value)) + geom_boxplot() +
    geom_hline(yintercept = 0, col = 'red', size = 1, linetype = 'dashed') +
    coord_flip() +
    labs(y = expression(paste(Delta, ' ', ln, ' ', p(y))), x = '') +
    theme_bw() +
    theme(axis.text  = element_text(size=25),
          axis.title = element_text(size=25,face="bold")) +
    scale_y_continuous(breaks = seq(-60, 0, 10))





