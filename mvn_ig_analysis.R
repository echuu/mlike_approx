

# path for lenovo
# setwd("C:/Users/ericc/mlike_approx")

# path for dell
setwd("C:/Users/chuu/mlike_approx")

library(mvtnorm)           # for draws from multivariate normal
library("numDeriv")        # for grad() function - numerical differentiation
library('MCMCpack')        # for rinvgamma() function

source("partition/partition.R")      # load partition extraction functions
source("mvn_ig_helper.R")  # load functions specific to this model


set.seed(1)

D_vec = c(3, 5, 7, 10)

n_cases = length(D_vec)
LIL_d   = numeric(n_cases)

N_approx = 1000 # number of approximations to compute

# store the approximations for each LIL in a column
approx_d = matrix(NA, N_approx, n_cases) # (N_approx x n_cases)

for (i in 1:n_cases) {
    
    D = D_vec[i]
    
    print(paste("d =", D))
    
    p = D - 1        # dimension of beta
    N = 50           # sample size
    
    mu_beta = rep(0, p)      # mu_beta
    V_beta = diag(1, p)      # V_beta
    a_0 = 2 / 2              # a
    b_0 = 1 / 2              # b 
    
    # beta    = c(5, -1)     # (p x 1) true coefficient vector
    beta    = sample(-10:10, p, replace = T)
    
    if (length(beta) != p) {
        print("dimension of beta must equal p")
    }
    
    sigmasq = 4            # (1 x 1) true variance
    
    I_N = diag(1, N)       # (N x N) identity matrix
    I_p = diag(1, p)       # (p x p) identity matrix
    
    X = matrix(rnorm(N * p), N, p) # (N x p) design matrix
    
    eps = t(rmvnorm(1, mean = rep(0, N), sigma = sigmasq * I_N)) # (N x 1)
    
    y = X %*% beta + eps # (N x 1) response vector
    
    
    # compute posterior parameters
    V_beta_inv = solve(V_beta)
    V_star_inv = t(X) %*% X + V_beta_inv
    
    V_star  =  solve(V_star_inv)                                  # (p x p)
    mu_star =  V_star %*% (t(X) %*% y + V_beta_inv %*% mu_beta)   # (p x 1)
    a_n     =  a_0 + N / 2                                        # (1 x 1)
    b_n     =  b_0 + 0.5 * (t(y) %*% y +                          # (1 x 1)
                                t(mu_beta) %*% V_beta_inv %*% mu_beta - 
                                t(mu_star) %*% V_star_inv %*% mu_star) %>%  c()
    
    # create prior, posterior objects
    prior = list(V_beta = V_beta, 
                 mu_beta = mu_beta, 
                 a_0 = a_0, 
                 b_0 = b_0,
                 y = y, X = X,
                 V_beta_inv = V_beta_inv)
    
    post  = list(V_star  =  V_star,
                 mu_star =  mu_star,
                 a_n     =  a_n,
                 b_n     =  b_n,
                 V_star_inv = V_star_inv)
    
    LIL_d[i] = lil(y, X, prior, post)
    
    # --------------------------------------------------------------------------
    
    
    def_approx = approx_lil(N_approx, prior, post, D) # (N_approx x 1)
    
    
    approx_d[,i] = def_approx
    
} # end loop


write.csv(approx_d, "sim_results.csv", row.names = F)

test_read = read.csv("sim_results.csv")


read.csv("sim_results.csv")

# rbind(LIL_d, colMeans(approx_d))

### plot true value as dashed line, plot approximation points

approx_df = data.frame(lil_hat = approx_d[,1], lil = LIL_d[1]) %>% 
    mutate(iter = 1:N_approx)

ggplot(approx_df, aes(iter, lil_hat)) + geom_point(col = 'blue') + 
    geom_hline(aes(yintercept = lil), 
               linetype = 'dashed', size = 1.5, color = "red")






