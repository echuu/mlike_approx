

# path for lenovo
setwd("C:/Users/ericc/mlike_approx")

# path for dell
# setwd("C:/Users/chuu/mlike_approx")

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
    
    
    # def_approx = approx_lil(N_approx, prior, post, D) # (N_approx x 1)
    
    
    # approx_d[,i] = def_approx
    
} # end loop

LIL_d

# different because theres more RNG in between the steps of computing the LIL
# for different values of D
# LIL_true = c(-110.9457, -116.7249, -140.0094, -129.2714)


write.csv(approx_d, "sim_results.csv", row.names = F)


# ------------------------------------------------------------------------------

# use these since they are the correct versions of the true LIL
LIL_true = c(-110.9457, -120.3410, -120.7169, -133.6093) # DELL OUTPUT 

approx_sim = read.csv("sim_results.csv")

approx_means = colMeans(approx_sim) %>% unlist %>% unname()

approx_med = apply(approx_sim, 2, median) %>%  unlist() %>% unname() 

approx_sds = apply(approx_sim, 2, sd) %>%  unlist() %>% unname() 

approx_means[1]

### plot true value as dashed line, plot approximation points

plot_approx = function(approx, approx_mean, approx_sd, LIL_true) {
    
    N_approx = length(approx)
    
    # (N_approx x 3) dataframe constructed so that we can use ggplot2
    approx_df = data.frame(lil_hat = approx, lil = LIL_true) %>% 
        mutate(iter = 1:N_approx)
    
    p = ggplot(approx_df, aes(iter, lil_hat)) + geom_point(col = 'blue') + 
        geom_hline(aes(yintercept = lil), 
                   linetype = 'dashed', size = 1.5, color = "red") +
        ggtitle(paste("LIL = ", LIL_true, 
                      ", Mean Approx = ", round(approx_mean, 3), 
                      ", SD Approx = ", round(approx_sd, 3),
                      sep = '')) + 
        theme(plot.title = element_text(size = 18))
    
    return(p)
}


plot_approx(approx_sim[,1], approx_means[1], approx_sds[1], LIL_true[1])

plot_approx(approx_sim[,2], approx_means[2], approx_sds[2], LIL_true[2])

plot_approx(approx_sim[,3], approx_means[3], approx_sds[3], LIL_true[3])

plot_approx(approx_sim[,4], approx_means[4], approx_sds[4], LIL_true[4])



# compare median vs mean for approximation results

approx_compare = rbind(LIL_true, approx_med) %>% rbind(approx_means)

# seems like the median does better for these approximations
colnames(approx_compare) = c("d = 3", "d = 5", "d = 7", "d = 10")

approx_compare






