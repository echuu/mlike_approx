
library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio

DELL_PATH = "C:/Users/chuu/mlike_approx"
# LEN_PATH  = "C:/Users/ericc/mlike_approx"
# path for lenovo
# setwd(LEN_PATH)

# path for dell
setwd(DELL_PATH)

source("partition/partition.R")
source("hybrid_approx.R")
source("mvn_ig/mvn_ig_helper.R") # load this LAST to overwrite def preprocess()



# STAN sampler settings --------------------------------------------------------

J         = 100          # number of MC samples per approximation
N_approx  = 10            # number of approximations
burn_in   = 2000         # number of burn in draws
n_chains  = 4            # number of markov chains to run
stan_seed = 123          # seed

J_iter = 1 / n_chains * N_approx * J + burn_in 


K_sims = 2               # num of simulations to run FOR EACH N in N_vec

# ------------------------------------------------------------------------------


D_vec = c(3, 5, 7, 10)
D_vec = c(3, 5)
LIL_d = vector("list", length = length(D_vec))    

set.seed(123)
for (d_i in 1:length(D_vec)) {
    
    D = D_vec[d_i]
    
    print(paste("Performing LIL approximation for D = ", D, sep = ''))
    
    # priors, true parameter values --------------------------------------------
    # D       = 3                # dimension of paramter
    p       = D - 1            # dimension of beta
    mu_beta = rep(0, p)        # prior mean for beta
    V_beta  = diag(1, p)       # scaled precision matrix for betaV_beta
    a_0     = 2 / 2            # shape param for sigmasq
    b_0     = 1 / 2            # scale param 
    beta    = sample(-10:10, p, replace = T)
    sigmasq = 4                # true variance (1 x 1) 
    
    I_p = diag(1, p)       # (p x p) identity matrix
    
    # --------------------------------------------------------------------------
    
    N_vec = c(50, 100, 150, 200, 300)
    N_vec = c(200) # for testing -- comment this line to perform ext. analysis
    
    
    LIL_N = numeric(length(N_vec))      # store the LIL for each N
    LIL_N_hat = numeric(length(N_vec))  # store LIL approximations for N
    
    # compute the true value of the LIL for each N 
    for (i in 1:length(N_vec)) {
        
        N = N_vec[i]
        
        I_N = diag(1, N)       # (N x N) identity matrix
        
        LIL_N_k = numeric(K_sims)     # store the true LIL for K_sims results
        LIL_N_k_hat = numeric(K_sims) # store the approx for K_sims results
        
        
        print(paste('iter = ', i, ' -- calculating LIL for N = ', N, 
                    ' (', K_sims, ' sims)', sep = ''))
        
        for (k in 1:K_sims) {
            
            ## simulate N data points + sample from posterior ------------------
            X = matrix(rnorm(N * p), N, p) # (N x p) design matrix
            
            eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))
            
            y = X %*% beta + eps # (N x 1) response vector
            # ------------------------------------------------------------------
            
            
            ## compute posterior parameters ------------------------------------
            V_beta_inv = solve(V_beta)
            V_star_inv = t(X) %*% X + V_beta_inv
            
            V_star  = solve(V_star_inv)                                # (p x p)
            mu_star = V_star %*% (t(X) %*% y + V_beta_inv %*% mu_beta) # (p x 1)
            a_n =  a_0 + N / 2 
            b_n =  c(b_0 + 0.5 * (t(y) %*% y + 
                                      t(mu_beta) %*% V_beta_inv %*% mu_beta - 
                                      t(mu_star) %*% V_star_inv %*% mu_star))
            
            # create prior, posterior objects
            prior = list(V_beta = V_beta, 
                         mu_beta = mu_beta, 
                         a_0 = a_0, 
                         b_0 = b_0,
                         y = y, X = X,
                         V_beta_inv = V_beta_inv)
            
            # store posterior parameters
            post  = list(V_star  =  V_star,
                         mu_star =  mu_star,
                         a_n     =  a_n,
                         b_n     =  b_n,
                         V_star_inv = V_star_inv)
            # ------------------------------------------------------------------
            
            ## compute true log marginal likelihood ----------------------------
            
            LIL_N_k[k] = lil(y, X, prior, post)
            
            # ------------------------------------------------------------------
            
            ## form the approximation
            post_dat = list(p = p,
                            a_n = a_n, b_n = b_n, 
                            mu_star = c(mu_star), V_star = V_star)
            
            mvnig_fit = stan(file   = 'mvn_ig/mvn_ig_sampler.stan', 
                             data    = post_dat,
                             iter    = J_iter,
                             warmup  = burn_in,
                             chains  = n_chains,
                             seed    = stan_seed,
                             refresh = 0) # should give us J * N_approx draws
            
            # use special preprocess b/c we call psi_true() 
            u_df = preprocess(mvnig_fit, D, post)
            
            # LIL_N_k_hat[k] = mean(approx_lil(N_approx, prior, D, u_df, J))
            LIL_N_k_hat[k] = mean(approx_lil(N_approx, D, u_df, J, prior))

        }
        
        LIL_N[i] = mean(LIL_N_k)
        LIL_N_hat[i] = mean(LIL_N_k_hat)

        
    } # end of main loop
    
    
    LIL_d[[d_i]] = rbind(LIL_N, LIL_N_hat)
}

# LIL_N_k
# LIL_N_k_hat

LIL_d

# output from old functions
# > LIL_d
# [[1]]
# [,1]
# LIL_N     -423.0917
# LIL_N_hat -421.1374
# 
# [[2]]
# [,1]
# LIL_N     -451.0669
# LIL_N_hat -444.3279





# ------------------------------------------------------------------------------


# verifying that LIL vs log N regression has -D/2 slope

# priors
D       = 6                # dimension of paramter
p       = D - 1            # dimension of beta
mu_beta = rep(0, p)        # prior mean for beta
V_beta  = diag(1, p)       # scaled precision matrix for betaV_beta
a_0     = 2 / 2            # shape param for sigmasq
b_0     = 1 / 2            # scale param 


# true model parameters  
beta    = sample(-10:10, p, replace = T)   # coefficient vector (p x 1)
sigmasq = 4                                # true variance (1 x 1) 



I_p = diag(1, p)       # (p x p) identity matrix

# values of N for which we will compute + approximate the LIL
N_vec = seq(50, 5000, 100)
LIL_N = numeric(length(N_vec)) # store the LIL for each of the grid values of N
K_sims = 200  # num of simulations to run FOR EACH N in N_vec

set.seed(1)
for (i in 1:length(N_vec)) {
    
    N = N_vec[i]
    
    I_N = diag(1, N)       # (N x N) identity matrix
    
    LIL_N_k = numeric(K_sims) # store the K_sims results
    
    for (k in 1:K_sims) {
        
        
        ## simulate N data points + sample from posterior ----------------------
        X = matrix(rnorm(N * p), N, p) # (N x p) design matrix
        
        eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))
        
        y = X %*% beta + eps # (N x 1) response vector
        # ----------------------------------------------------------------------
        
        
        ## compute posterior parameters ----------------------------------------
        V_beta_inv = solve(V_beta)
        V_star_inv = t(X) %*% X + V_beta_inv
        
        V_star  =  solve(V_star_inv)                                  # (p x p)
        mu_star =  V_star %*% (t(X) %*% y + V_beta_inv %*% mu_beta)   # (p x 1)
        a_n =  a_0 + N / 2 
        b_n =  b_0 + 0.5 * (t(y) %*% y + 
                                t(mu_beta) %*% V_beta_inv %*% mu_beta - 
                                t(mu_star) %*% V_star_inv %*% mu_star) %>%  c()
        
        # create prior, posterior objects
        prior = list(V_beta = V_beta, 
                     mu_beta = mu_beta, 
                     a_0 = a_0, 
                     b_0 = b_0,
                     y = y, X = X,
                     V_beta_inv = V_beta_inv)
        
        # store posterior parameters
        post  = list(V_star  =  V_star,
                     mu_star =  mu_star,
                     a_n     =  a_n,
                     b_n     =  b_n,
                     V_star_inv = V_star_inv)
        # ----------------------------------------------------------------------
        
        ## compute true log marginal likelihood --------------------------------
        
        # beta_mle = solve(t(X) %*% X) %*% t(X) %*% y
        ybar = X %*% mu_star
        sigmasq_mle = 1 / N * sum((y - ybar)^2)
        
        LIL_N_k[k] = lil(y, X, prior, post) - 
            sum(dnorm(y, ybar, sqrt(sigmasq_mle), log = T))
        
        # a_0 * log(b_0) + lgamma(a_n) + 0.5 * log_det(V_star) - 
        #    0.5 * log_det(V_star) - N / 2 * log(2 * pi) - lgamma(a_0) - 
        #    a_n * log(b_n)
        
        
        
    }
    
    LIL_N[i] = mean(LIL_N_k)
    
}

LIL_N



LIL_df = data.frame(LIL_N = LIL_N, log_N = log(N_vec))
plot(LIL_N ~ log_N, LIL_df)

lm(LIL_N ~ log_N, LIL_df)












