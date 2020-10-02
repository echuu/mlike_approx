

setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions
setwd("C:/Users/ericc/mlike_approx/covariance")
source("HIW_helper.R")      # covariance related helper functions
source("HIW_graph.R")


testG  = matrix(c(1,1,0,0,0,
                  1,1,1,1,1,
                  0,1,1,1,1,
                  0,1,1,1,1,
                  0,1,1,1,1), 5, 5)

D = nrow(testG)
b = 3          # prior degrees of freedom
V = diag(1, D) # prior scale matrix 

D_0 = 0.5 * D * (D + 1) # num entries on diagonal and upper diagonal
J = 30



log_density = function(u, data) {
    -psi(u, data)
}

hme_exp_term = function(u) {
    Lt = matrix(0, D, D)     # (D x D) lower triangular matrix
    Lt_vec_0 = numeric(D_0)  # (D_0 x 1) vector to fill upper triangular, Lt
    Lt_vec_0[edgeInd] = u
    Lt[upper.tri(Lt, diag = T)] = Lt_vec_0   # populate lower triangular terms
    
    - N * log_det(Lt) + 0.5 * matrix.trace(t(Lt) %*% Lt %*% S)
    
}





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
N = 100
Y = matrix(0, N, D)
for (i in 1:N) {
    Y[i, ] = t(t(chol(Sigma_G)) %*% rnorm(D, 0, 1)) # (500 x D)
}

S = t(Y) %*% Y

params = list(N = N, D = D, D_0 = D_0, testG = testG, edgeInd = edgeInd,
              upperInd = upperInd, S = S, V = V, b = b)

postIW = sampleHIW(J, D_u, D_0, testG, b, N, V, S, edgeInd)  
post_samps = postIW$post_samps                 # (J x D_u)

u_df = preprocess(post_samps, D_u, params)     # J x (D_u + 1)

(LIL = logmarginal(Y, testG, b, V, S))

# (gnorm_approx = - 0.5 * D * N * log(2 * pi) + 
#         gnorm(testG, b + N, V + S, iter = 100) - gnorm(testG, b, V, iter = 100))

hybrid = hybrid_ml(D_u, u_df, J, params)
hybrid$zhat

u_samp = as.matrix(post_samps)
colnames(u_samp) = names(u_df)[1:D_u]
# prepare bridge_sampler input()
lb = rep(-Inf, D_u)
ub = rep(Inf, D_u)
names(lb) <- names(ub) <- colnames(u_samp)

bridge_result = bridgesampling::bridge_sampler(samples = u_samp, 
                                               log_posterior = log_density,
                                               data = params, 
                                               lb = lb, ub = ub, 
                                               silent = TRUE)
bridge_result$logml


hme_approx = function(u_df, params, J, D, N) {
    
    N   = params$N
    D   = params$D
    D_0 = params$D_0
    S   = params$S
    
    tmp_J = apply(u_df[,1:D_u], 1, hme_exp_term) # term in the exponential
    log(J) - 0.5 * N * D * log(2 * pi) - log_sum_exp(tmp_J)
    
}

hme_approx(u_df, params, J, D, N)


aB = 100 # number of replications
hyb = numeric(B)
bridge = numeric(B)
set.seed(123)

for (i in 1:B) {
    
    postIW = sampleHIW(J, D_u, D_0, testG, b, N, V, S, edgeInd)  
    post_samps = postIW$post_samps                 # (J x D_u)
    
    u_df = preprocess(post_samps, D_u, params)     # J x (D_u + 1)
    
    hybrid = hybrid_ml(D_u, u_df, J, param_list)

    if (any(is.na(hybrid))) {print(paste("error in iteration", i)); next;}
    
    hyb[i] = hybrid$zhat
    
    u_samp = as.matrix(post_samps)
    colnames(u_samp) = names(u_df)[1:D_u]
    lb = rep(-Inf, D_u)
    ub = rep(Inf, D_u)
    names(lb) <- names(ub) <- colnames(u_samp)
    bridge_result = bridgesampling::bridge_sampler(samples = u_samp, 
                                                   log_posterior = log_density,
                                                   data = params, 
                                                   lb = lb, ub = ub, 
                                                   silent = TRUE)
    bridge[i] = bridge_result$logml
    
    
    avg_hyb = mean(hyb[hyb!=0])
    avg_bridge = mean(bridge[bridge!=0])
    
    print(paste("iter ", i, ': ',
                "hybrid = ", round(avg_hyb, 3), 
                '; ', "mae = ", round(mean(abs(LIL - hyb[hyb!=0])), 4),
                ' // ',
                "bridge = ", round(avg_bridge, 3), '; ', 
                "mae = ", round(mean(abs(LIL - bridge[bridge!=0])), 4),
                sep = '')) 
}
approx = data.frame(LIL, hyb = hyb[hyb!=0], bridge = bridge[bridge!=0])
error = data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           mae = colMeans(abs(LIL - approx)),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)


saveRDS(list(J = J, D = D, D_u = D_u, N = N, approx_df = approx, error = error), 
        file = 'covarIW_d5_j5000.RData')
mvnig_d20 = readRDS('covarIW_d5.RData')







