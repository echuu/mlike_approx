

setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions
setwd("C:/Users/ericc/mlike_approx/covariance")

source("HIW_graph.R")
source("misc_graph_funcs.R")
source("HIW_helper.R")


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

testG = matrix(c(1,1,1,0,0,
                 1,1,1,0,0,
                 1,1,1,1,1,
                 0,0,1,1,1,
                 0,0,1,1,1), 5,5)
D = nrow(testG)

b = 3          # prior degrees of freedom
V = diag(1, D) # prior scale matrix 

D_0 = 0.5 * D * (D + 1) # num entries on diagonal and upper diagonal
J = 2000


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

postIW = sampleHIW(J, D_u, D_0, testG, b, N, V, S, edgeInd)  
post_samps = postIW$post_samps                 # (J x D_u)

logmarginal = function(Y, G, b, D, S){
    n = nrow(Y); p = ncol(Y); # S = t(Y)%*%Y
    logmarginal = -0.5*n*p*log(2*pi) + logHIWnorm(G, b, D) - 
        logHIWnorm(G, b+n, D+S)
    return(logmarginal)
}

## testing ---------------------------------------------------------------------


u_df = preprocess(post_samps, D_u, params)     # J x (D_u + 1)
hml_approx = hml_const(1, D_u, u_df, J, params)
gnorm_approx = - 0.5 * D * N * log(2 * pi) + gnorm(testG, b + N, V + S, iter = 100) - 
    gnorm(testG, b, V, iter = 100)

logmarginal(Y, testG, b, V, S)
gnorm_approx
hml_approx$const_vec

## bridge

log_density = function(u, data) {
    -psi(u, data)
}

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


## bridge ----------------------------------------------------------------------

B = 100
J = 1000
n_samps = 10


hyb     = numeric(B) # store hybrid estiator
hyb_fs  = numeric(B) # store hybrid fs estiator
hyb_ss  = numeric(B) # store harmonic mean estiator
hyb_ts  = numeric(B) # store harmonic mean estiator
bridge  = numeric(B)

resample_partition = function(fit_partition, D, n_resamp) {
    out <- tryCatch(
        {
            fit_resid(fit_partition, D, n_samps, params)
            
        },
        error=function(cond) {
            message(paste("resample computation error"))
            message(cond)
            return(NA)
        },
        warning=function(cond) {
            message("Here's the original warning message:")
            message(cond)
            return(NULL)
        },
        finally={
            
        }
    )    
    return(out)
}
hybrid <- function(u_df) {
    out <- tryCatch(
        {
            
            hml_const(1, ncol(u_df) - 1, u_df, J, params)
            # The return value of `readLines()` is the actual value 
            # that will be returned in case there is no condition 
            # (e.g. warning or error). 
            # You don't need to state the return value via `return()` as code 
            # in the "try" part is not wrapped insided a function (unlike that
            # for the condition handlers for warnings and error below)
        },
        error=function(cond) {
            message(paste("hybrid computation error"))
            message(cond)
            return(NA)
        },
        warning=function(cond) {
            message("Here's the original warning message:")
            message(cond)
            return(NULL)
        },
        finally={
            
        }
    )    
    return(out)
}

set.seed(1)
for (b_i in 1:B) {
    
    postIW = sampleHIW(J, D_u, D_0, testG, b, N, V, S, edgeInd)  
    post_samps = postIW$post_samps                 # (J x D_u)
    
    u_df = preprocess(post_samps, D_u, params)     # J x (D_u + 1)
    hml_approx = hml_conshit(1, D_u, u_df, J, params)
    hml_approx = hybrid(u_df)
    if (any(is.na(hml_approx))) {print(paste("error in iteration", b_i)); next;}
    
    og_part = hml_approx$param_out %>% dplyr::select(-c(psi_choice, logQ_cstar))
    ss_part = resample_partition(og_part, D_u, n_samps)
    if (any(is.na(ss_part))) {print(paste("error in iteration", b_i)); next;}
    
    ts_part = resample_partition(ss_part, D_u, n_samps / 2)
    if (any(is.na(ts_part))) {print(paste("error in iteration", b_i)); next;}
    
    hyb_fs[b_i] = hml_approx$const_vec
    hyb_ss[b_i] = log_sum_exp(unlist(compute_expterms(ss_part, D_u)))
    hyb_ts[b_i] = log_sum_exp(unlist(compute_expterms(ts_part, D_u)))
    
    hyb[b_i] = (hyb_fs[b_i] + hyb_ss[b_i] + hyb_ts[b_i]) / 3
    
    u_samp = as.matrix(post_samps)
    colnames(u_samp) = names(u_df)[1:D_u]

    bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
                                                   log_posterior = log_density,
                                                   data = params,
                                                   lb = lb, ub = ub,
                                                   silent = TRUE)
    bridge[b_i] = bridge_result$logml
    
    
    print(paste("iter: ", b_i, '; ',
                'bridge: ', round(mean(bridge[1:b_i]), 3), '; ',
                'hybrid: ', round(mean(hyb[hyb!=0]), 4),
                sep = '')) 
    
    # if (b_i %% 1 == 0) { 
    #     print(paste("iter: ", b_i, '; ',
    #                 'bridge: ', round(mean(bridge[1:b_i]), 3), '; ',
    #                 'hybrid: ', round(mean(hyb[1:b_i]), 4),
    #                 sep = '')) 
    # }
    
}


LIL = logmarginal(Y, testG, b, V, S)
approx = data.frame(bridge = bridge, 
                    hyb = hyb[hyb!=0], 
                    hyb_fs = hyb_fs[hyb_fs!= 0],
                    hyb_ss = hyb_ss[hyb_ss!= 0], 
                    hyb_ts = hyb_ts[hyb_ts!=0], LIL)

# approx = data.frame(hyb[1:B], came[1:B], came_0)
data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           ae = colMeans(LIL - approx),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)





