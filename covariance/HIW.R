

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
J = 2000


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

# logical vector determining existence of edges between vertices
edgeInd = testG[upper.tri(testG, diag = TRUE)] %>% as.logical
upperInd = testG[upper.tri(testG)] %>% as.logical
D_u = sum(edgeInd)

# Specify true value of Sigma
set.seed(1)
true_params = HIWsim(testG, b, V)
Sigma_G = true_params$Sigma
Omega_G = true_params$Omega

# chol(Omega_G)

# Generate data Y based on Sigma_G
N = 200
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

logmarginal = function(Y, G, b, D, S){
    n = nrow(Y); p = ncol(Y); # S = t(Y)%*%Y
    logmarginal = -0.5*n*p*log(2*pi) + logHIWnorm(G, b, D) - 
        logHIWnorm(G, b+n, D+S)
    return(logmarginal)
}


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
hml_approx$const_vec

hml_approx$param_out

abs(logmarginal(Y, testG, b, V, S) - hml_approx$const_vec)

abs(logmarginal(Y, testG, b, V, S) - gnorm_approx)

hml_approx$param_out %>% 
    dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs)


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












