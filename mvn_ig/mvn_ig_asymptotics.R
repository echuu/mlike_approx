
library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio

# DELL_PATH = "C:/Users/chuu/mlike_approx"
LEN_PATH  = "C:/Users/ericc/mlike_approx"
# path for lenovo
setwd(LEN_PATH)

# path for dell
# setwd(DELL_PATH)

source("partition/partition.R")
source("extractPartition.R")
source("hybrid_approx.R")
source("mvn_ig/mvn_ig_helper.R") # load this LAST to overwrite def preprocess()



# STAN sampler settings --------------------------------------------------------

J         = 500          # number of MC samples per approximation
N_approx  = 1            # number of approximations
burn_in   = 2000         # number of burn in draws
n_chains  = 4            # number of markov chains to run
stan_seed = 123          # seed

J_iter = 1 / n_chains * N_approx * J + burn_in 


K_sims = 1               # num of simulations to run FOR EACH N in N_vec

# ------------------------------------------------------------------------------


# D_vec = c(3, 5, 7, 10)
D_vec = c(3)
LIL_d = vector("list", length = length(D_vec))    


LIL_const  = matrix(NA, N_approx, length(N_vec))
LIL_hybrid = matrix(NA, N_approx, length(N_vec))
LIL_taylor  = matrix(NA, N_approx, length(N_vec))


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
    
    # N_vec = c(50, 60, 70, 100, 110, 125, 150, 200, 225, 250, 300)
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
            
            # compute MLE estimates for mean, variance of the regression model
            ybar = X %*% mu_star
            sigmasq_mle = 1 / N * sum((y - ybar)^2)
            
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
            
            # subtract log max likelihood to stabilize approximation
            LIL_N_k[k] = lil(y, X, prior, post) #- 
                #sum(dnorm(y, ybar, sqrt(sigmasq_mle), log = T))
            
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
                             refresh = 0) # should give us J * N_approx draws
            
            # use special preprocess b/c we call psi_true() 
            u_df = preprocess(mvnig_fit, D, post)
            
            u_df %>% head
            
            # subtract log max likelihood to stabilize approximation
            LIL_N_k_hat[k] = mean(approx_lil(N_approx, D, u_df, J, prior)) #- 
                #sum(dnorm(y, ybar, sqrt(sigmasq_mle), log = T))
            
            hml_approx = hml(N_approx, D, u_df, J, prior) 
            
            LIL_const[,k]  = hml_approx$const_vec
            LIL_taylor[,k] = hml_approx$taylor_vec
            LIL_hybrid[,k] = hml_approx$hybrid_vec
            
            # comment out later 
            hml_approx$const_vec   # -429.8105
            hml_approx$taylor_vec  # -427.3445
            hml_approx$hybrid_vec  # -427.3445
            
            hml_approx$n_taylor
            hml_approx$n_const
            

        } # end of K_sims loop
        
        LIL_N[i] = mean(LIL_N_k)
        LIL_N_hat[i] = mean(LIL_N_k_hat)
        
        print(rbind(LIL_N[i], LIL_N_hat[i]))

        
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
D       = 3                # dimension of paramter
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
N_vec = c(50, 60, 70, 100, 110, 125, 150, 200, 225, 250, 300)
# N_vec = c(200)
LIL_N = numeric(length(N_vec)) # store the LIL for each of the grid values of N
K_sims = 100  # num of simulations to run FOR EACH N in N_vec

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


# ------------------------------------------------------------------------------
library(reshape2)

N_vec = c(50, 60, 70, 100, 110, 125, 150, 200, 225, 250, 300)


# D = 3
LIL_N = c(-116.3446, -137.6999, -159.6706, -222.8907, -244.2451, -275.9732,
          -329.4240, -433.5548, -488.2828, -542.3228, -648.5599)

LIL_N_hat = c(-103.9668, -125.1452, -147.3270, -210.9323, -232.3533, -264.2786,
              -317.4581, -421.9521, -476.5969, -530.7302, -637.8043)


lil_df = data.frame(LIL_N = LIL_N, LIL_N_hat = LIL_N_hat, log_N = log(N_vec))

lil_df_long = melt(lil_df, id.vars = "log_N")

ggplot(lil_df_long, aes(log_N, value, col = variable)) + geom_point()

# D = 5
LIL_N = c(-131.1138, -153.8555, -175.8266, -240.4139, -261.4032, -294.8986, 
          -348.9524, -455.3916, -508.6343, -560.8775, -667.1969)

LIL_N_hat = c(-111.8813, -135.1093, -157.1889, -222.3254, -243.5304, -276.9260, 
              -331.3026, -437.6792, -491.3586, -543.9151, -651.0929)


# D = 7
LIL_N = c(-135.8489, -158.0649, -181.5670, -247.1824, -267.3347, -300.0477,
          -354.3124, -460.4378, -516.1904, -566.7185, -674.6706)

LIL_N_hat = c(-110.1174, -132.1217, -156.2267, -222.1395, -242.8230, -275.8084,
              -330.0965, -436.1873, -492.7829, -543.5327, -652.6138)


# D = 10
LIL_N = c(-147.2810, -170.2424, -193.1394, -260.4092, -282.9458, -315.8182,
          -369.3659, -478.8867, -532.4860, -586.4090, -693.5570)

LIL_N_hat = c(-109.8516, -133.9974, -157.3436, -225.7665, -249.0499, -281.6845,
              -335.8194, -445.8283, -499.9098, -554.2056, -662.9794)



