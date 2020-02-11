

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



# sampler settings -------------------------------------------------------------

N_approx  = 1            # number of approximations
burn_in   = 2000         # number of burn in draws
n_chains  = 4            # number of markov chains to run
stan_seed = 123          # seed
K_sims    = 20           # num of simulations to run FOR EACH N in N_vec


# ------------------------------------------------------------------------------

J_samps = c(40, 100, 160, 200, 240, 300, 500, 1000, 2000, 5000)
D_vec   = c(3, 5, 8, 10, 15, 20, 25, 30, 40, 50)

# J_samps = c(40, 100)
# D_vec = 3

N = 200 # fix this value for now


LIL_const_list  = vector("list", length = length(D_vec))    
LIL_taylor_list = vector("list", length = length(D_vec))    
LIL_hybrid_list = vector("list", length = length(D_vec))    
LIL_N_list      = vector("list", length = length(D_vec))


set.seed(123)

for (d in 1:length(D_vec)) {
    
    
    D = D_vec[d]
    p       = D - 1            # dimension of beta
    mu_beta = rep(0, p)        # prior mean for beta
    V_beta  = diag(1, p)       # scaled precision matrix for betaV_beta
    a_0     = 2 / 2            # shape param for sigmasq
    b_0     = 1 / 2            # scale param 
    beta    = sample(-10:10, p, replace = T)
    sigmasq = 4                # true variance (1 x 1) 
    I_p = diag(1, p)           # (p x p) identity matrix
    

    LIL_const  = matrix(NA, K_sims, length(J_samps))
    LIL_hybrid = matrix(NA, K_sims, length(J_samps))
    LIL_taylor = matrix(NA, K_sims, length(J_samps))
    LIL_N_k    = matrix(NA, K_sims, length(J_samps))
    
    for (i in 1:length(J_samps)) {
        
        J = J_samps[i] # number of MCMC samples to use in approximation
        
        J_iter = 1 / n_chains * N_approx * J + burn_in 
        
        print(paste('D = ', D, ': iter = ', i, '/', length(J_samps), ' -- J = ', 
                    J, ' MCMC samples (', K_sims, ' sims)', sep = ''))
        
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
            LIL_N_k[k,i] = lil(y, X, prior, post)
            # ------------------------------------------------------------------
            
            # ALGORITHM BEGINS HERE
            
            ## sample from posterior
            post_dat = list(p = p,
                            a_n = a_n, b_n = b_n, 
                            mu_star = c(mu_star), V_star = V_star)
            
            mvnig_fit = stan(file   = 'mvn_ig/mvn_ig_sampler.stan', 
                             data    = post_dat,
                             iter    = J_iter,
                             warmup  = burn_in,
                             chains  = n_chains,
                             refresh = 0) # should give us J * N_approx draws
            
            ## evaluate psi(U)
            u_df = preprocess(mvnig_fit, D, post)
            
            ## compute approximation (3 different ways)
            
            hml_approx = hml(N_approx, D, u_df, J, prior) 
            
            # Note: only 1 approximation per method per simulation computed
            LIL_const[k,i]  = hml_approx$const_vec    
            LIL_taylor[k,i] = hml_approx$taylor_vec
            LIL_hybrid[k,i] = hml_approx$hybrid_vec
            
        }
        
    } # end of J_samps loop
    
    LIL_const_list[[d]]  = LIL_const
    LIL_taylor_list[[d]] = LIL_taylor
    LIL_hybrid_list[[d]] = LIL_hybrid
    LIL_N_list[[d]]      = LIL_N_k
    
    print(do.call(rbind, list(colMeans(LIL_N_k),  
                              colMeans(LIL_const),  
                              colMeans(LIL_taylor),  
                              colMeans(LIL_hybrid))))

} # end of D_vec loop


LIL_N_k
LIL_const
LIL_taylor
LIL_hybrid


LIL_const_list
LIL_taylor_list
LIL_hybrid_list
LIL_N_list

require(rlist)
list.save(LIL_const_list,  'const.rds')
list.save(LIL_taylor_list, 'taylor.rds')
list.save(LIL_hybrid_list, 'hybrid.rds')
list.save(LIL_N_list,      'true.rds')

test_read = list.load('taylor.rds')


colMeans(LIL_N_k)
colMeans(LIL_const)
colMeans(LIL_taylor)
colMeans(LIL_hybrid)

do.call(rbind, list(colMeans(LIL_N_list[[1]]),  colMeans(LIL_const_list[[1]]),  
                    colMeans(LIL_taylor_list[[1]]),  
                    colMeans(LIL_hybrid_list[[1]])))


# ------------------------------------------------------------------------------



for (i in 1:length(D_vec)) {
    
    
    J_df = do.call(rbind, list(colMeans(LIL_N_list[[i]]),  
                               colMeans(LIL_const_list[[i]]),  
                               colMeans(LIL_taylor_list[[i]]),  
                               colMeans(LIL_hybrid_list[[i]]))) %>% data.frame
    
    names(J_df) = paste(rep('J = '), J_samps, sep = '')
    row.names(J_df) = c("LIL", "const", "taylor", "hybrid")
    
}



















