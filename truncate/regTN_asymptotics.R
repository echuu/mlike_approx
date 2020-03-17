
library(TruncatedNormal)
library(tmg)
library(mvtnorm)

LEN_PATH  = "C:/Users/ericc/mlike_approx"
setwd(LEN_PATH)

source("partition/partition.R")
source("extractPartition.R")
source("hybrid_approx.R")
source("truncate/regTN_helper.R")






# ------------------------------------------------------------------------------

set.seed(1)
K_sims    = 100                              # num of sims to run for each N
N_vec_log = seq(5, 10, by = 0.25)            # sample size grid unif in log
N_vec     = floor(exp(N_vec_log)) %>% unique # sample size to generate data


LIL_N = numeric(length(N_vec))      # store the LIL for each N

for (i in 1:length(N_vec)) {
    
    N = N_vec[i]
    print(paste("N = ", N, sep = ''))
    
    
    LIL_N_k = numeric(K_sims)     # store the true LIL for K_sims
    LIL_N_k_hat = numeric(K_sims) # store the approx for K_sims results
    
    for (k in 1:K_sims) {
        
        # generate data --------------------------------------------------------
        X   = matrix(rnorm(N * D), N, D)                # (N x D) design matrix
        eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))    # (N x 1) errors vector
        y   = X %*% beta + eps                          # (N x 1) response vector
        
        # compute posterior parameters -------------------------------------------------
        Q_beta     =  1 / sigmasq * (t(X) %*% X + tau * diag(1, D))
        Q_beta_inv =  solve(Q_beta)
        b          =  1 / sigmasq * t(X) %*% y
        mu_beta    =  Q_beta_inv %*% b
        
        # create prior, post objects to be passed into the hml algorithm 
        prior = list(y = y, X = X, sigmasq = sigmasq, tau = tau, N = N, D = D)
        post = list(Q_beta = Q_beta, Q_beta_inv = Q_beta_inv, mu_beta = mu_beta, b = b)
        # ----------------------------------------------------------------------
        
        # compute true LIL first
        # TODO: without forming approximation, plot the true LIL to verify
        # -D/2 slope
        lil_0 = -0.5 * N * log(2 * pi) - 0.5 * (N + D) * log(sigmasq) + 
            0.5 * D * log(tau) - 0.5 * log_det(Q_beta) - 
            1 / (2 * sigmasq) * sum(y^2) + 0.5 * sum(b * mu_beta)
        
        LIL_N_k[k] = lil_0 +
            log(TruncatedNormal::pmvnorm(mu_beta, Q_beta_inv,
                                         lb = rep(0, D), ub = rep(Inf, D))[1]) -
            sum(dnorm(y, X %*% mu_beta, sqrt(sigmasq), log = T))
        
        # LIL_N_k[k] = lil_0 - 
        #     sum(dnorm(y, X %*% mu_beta, sqrt(sigmasq), log = T))

        
        # sample from posterior, u
        
        
        # evaluate psi(u)
        
        
        # compute approximation
        
        # subtract the loglikelihood
        
        
    } # end of k-sims
    
    
    LIL_N[i] = mean(LIL_N_k)
    
    
    # print(paste("approximating LIL for N = ", N, " -- LIL = ", 
    #             round(mean(LIL_N_k[i, k]), 2), sep = ''))
    
    
}



LIL_df = data.frame(LIL_N = LIL_N, log_N = log(N_vec))

LIL_df = melt(LIL_df, id.vars = "log_N")


library(reshape2)
library(ggpmisc)
formula1 = y ~ x

ggplot(LIL_df, aes(x = log_N, y = LIL_N)) + geom_point(size = 1.5) + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    labs(x = "log(n)", y = "log(ML)", 
         title = "Approximate ML for Logistic Regression, D = 3") + 
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16) + 
    theme(legend.position = "none")


ggplot(LIL_df, aes(x = log_N, y = value, 
                        color = as.factor(variable))) + geom_point(size = 1.3) + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    labs(x = "log(n)", y = "log(Z)", 
         title = "True (Red), Hybrid (Blue), Taylor (Green)") + 
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16) + 
    theme(legend.position = "none")







