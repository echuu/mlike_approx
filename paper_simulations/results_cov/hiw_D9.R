

# setwd("C:/Users/ericc/mlike_approx/algo")
# source("setup.R")           # setup global environment, load in algo functions
# setwd("C:/Users/ericc/mlike_approx/covariance")
# source("HIW_helper.R")      # covariance related helper functions
# source("HIW_graph.R")

setwd("/home/grad/ericchuu/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions
setwd("/home/grad/ericchuu/mlike_approx/covariance")
source("HIW_helper.R")      # covariance related helper functions
source("HIW_graph.R")

D = 9

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

# testG  = matrix(c(1,1,0,0,0,
#                   1,1,1,1,1,
#                   0,1,1,1,1,
#                   0,1,1,1,1,
#                   0,1,1,1,1), 5, 5)

D = nrow(testG)
b = 3          # prior degrees of freedom
V = diag(1, D) # prior scale matrix 

D_0 = 0.5 * D * (D + 1) # num entries on diagonal and upper diagonal
J = 5000


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

u_df = preprocess(post_samps, D_u, params)     # J x (D_u + 1)

logmarginal(Y, testG, b, V, S)

hml_approx = hml_const(1, D_u, u_df, J, params)
hml_approx$const_vec

B = 10 # number of replications
hyb_wt1  = numeric(B)  # store harmonic mean estiator
hyb_wt2  = numeric(B)  # store harmonic mean estiator
hyb_avg  = numeric(B)  # store harmonic mean estiator

for (b_i in 1:B) {
    
    postIW = sampleHIW(J, D_u, D_0, testG, b, N, V, S, edgeInd)  
    post_samps = postIW$post_samps                 # (J x D_u)
    
    u_df = preprocess(post_samps, D_u, params)     # J x (D_u + 1)
    
    hybrid = logml(D_u, u_df, J, params)
    if (any(is.na(hybrid))) {print(paste("error in iteration", b_i)); next;}
    
    hybrid$all_approx
    hybrid$wt_approx1
    hybrid$wt_approx2
    hybrid$wt_approx3
    
    hyb_wt1[b_i] = hybrid$wt_approx1
    hyb_wt2[b_i] = hybrid$wt_approx2
    hyb_wt3[b_i] = hybrid$wt_approx3
 
    print(paste("iter ", b_i, ': ',
                "hybrid_wt1 = ", round(mean(hyb_wt1[hyb_wt1!=0]), 3), '; ',
                "hybrid_wt2 = ", round(mean(hyb_wt2[hyb_wt2!=0]), 3), '; ',
                "hybrid_wt3 = ", round(mean(hyb_wt3[hyb_wt3!=0]), 3), '; ',
                sep = '')) 
       
}

LIL = logmarginal(Y, testG, b, V, S)
approx = data.frame(LIL, 
                    hyb_wt1 = hyb_wt1[hyb_wt1!=0], 
                    hyb_wt2 = hyb_wt2[hyb_wt2!=0], 
                    hyb_wt3 = hyb_wt3[hyb_wt3!=0])
data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           ae = colMeans(LIL - approx),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)









