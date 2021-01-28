

library(TruncatedNormal)
library(tmg)
library(mvtnorm)

library(Rcpp)
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")     
source("C:/Users/ericc/mlike_approx/truncate/regTN_helper.R")

sourceCpp("C:/Users/ericc/mlike_approx/speedup/trunc_psi.cpp")

l1_norm = function(u, u_0) {
    sum(abs(u - u_0))
}


set.seed(123)
D = 20
N = 100
I_D = diag(1, D)

n_samps = 10
J       = 500
B       = 100 # number of replications



# prior mean
mu_0 = rep(0, D)
tau     = 1 / 4          # precision: inverse of variance
sigmasq = 4              # true variance (1 x 1) 

# true value of beta
set.seed(1)
beta = sample(0:10, D, replace = T)
beta = runif(D)
beta = c(runif(D-1, 0, 1), 0)

# generate the regression data -------------------------------------------------

X   = matrix(rnorm(N * D), N, D)                # (N x D) design matrix
eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))    # (N x 1) errors vector
y   = X %*% beta + eps                          # (N x 1) response vector



# compute posterior parameters -------------------------------------------------
Q_beta     =  1 / sigmasq * (t(X) %*% X + tau * diag(1, D))
Q_beta_inv =  solve(Q_beta)
b          =  1 / sigmasq * t(X) %*% y
mu_beta    =  Q_beta_inv %*% b


# create prior, post objects to be passed into the hml algorithm 
prior = list(y = y, X = X, sigmasq = sigmasq, tau = tau, N = N, D = D,
             Q_beta = Q_beta, b = b)
post = list(Q_beta = Q_beta, Q_beta_inv = Q_beta_inv, mu_beta = mu_beta, b = b)



#### compute true log marginal likelihood --------------------------------------

lil_0 = -0.5 * N * log(2 * pi) - 0.5 * (N + D) * log(sigmasq) + 
    0.5 * D * log(tau) - 0.5 * log_det(Q_beta) - 
    1 / (2 * sigmasq) * sum(y^2) + 0.5 * sum(b * mu_beta)

(LIL  = lil_0 + D * log(2) + 
        log(TruncatedNormal::pmvnorm(mu_beta, Q_beta_inv, 
                                     lb = rep(0, D), ub = rep(Inf, D))[1]))

samples = rtmvnorm(J, c(mu_beta), Q_beta_inv, rep(0, D), rep(Inf, D))
u_df = preprocess(data.frame(samples), D, prior)

# write.csv(u_df, "C:/Users/ericc/u_df.csv", row.names = FALSE)

hybrid = hybrid_ml(D, u_df, J, prior)
hybrid$zhat


lb = rep(0, D)
ub = rep(Inf, D)
colnames(samples) = names(u_df)[1:D]
names(lb) <- names(ub) <- colnames(samples)

log_density = function(u, data) {
    -psi(u, data)
}

library(bridgesampling)
bridge_result <- bridge_sampler(samples = samples, log_posterior = log_density,
                                data = prior, lb = lb, ub = ub, silent = TRUE)
bridge_result$logml

LIL


ep = list(H_k = Q_beta, H_k_inv = Q_beta_inv, m_k = mu_beta, lambda = lambda)

hyb(u_df, ep)

hyb_numer(u_df, psi = psi_trunc)


set.seed(1)
D_MAX = 100
D_vec = seq(10, D_MAX, 2)
N_REP = 50
KAPPA = 3 # (# MCMC samples) / D = J / D

logml_truth = numeric(length(D_vec))
approx = matrix(0, N_REP, length(D_vec))
bridge = matrix(0, N_REP, length(D_vec))

i = 1
for (D in D_vec) {
    
    J = D * KAPPA
    
    ## priors ------------------------------------------------------------------
    mu_0 = rep(0, D)      # prior mean for beta
    tau  = 1 / 4          # precision: inverse of variance
    sigmasq = 4           # true variance (1 x 1) 
    
    ## true beta ---------------------------------------------------------------
    beta = sample(0:10, D, replace = T)
    beta = runif(D)
    beta = c(runif(D-1, 0, 1), 0)
    
    # generate the regression data ---------------------------------------------
    X   = matrix(rnorm(N * D), N, D)                # (N x D) design matrix
    eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))    # (N x 1) errors vector
    y   = X %*% beta + eps                          # (N x 1) response vector
    
    # compute posterior parameters ---------------------------------------------
    Q_beta     =  1 / sigmasq * (t(X) %*% X + tau * diag(1, D))
    Q_beta_inv =  solve(Q_beta)
    b          =  1 / sigmasq * t(X) %*% y
    mu_beta    =  Q_beta_inv %*% b
    
    
    # create prior, post objects to be passed into the hml algorithm 
    prior = list(y = y, X = X, sigmasq = sigmasq, tau = tau, N = N, D = D,
                 Q_beta = Q_beta, b = b)
    post = list(Q_beta = Q_beta, Q_beta_inv = Q_beta_inv, 
                mu_beta = mu_beta, b = b)
    
    #### compute true log marginal likelihood ----------------------------------
    lil_0 = -0.5 * N * log(2 * pi) - 0.5 * (N + D) * log(sigmasq) + 
        0.5 * D * log(tau) - 0.5 * log_det(Q_beta) - 
        1 / (2 * sigmasq) * sum(y^2) + 0.5 * sum(b * mu_beta)
    
    LIL  = lil_0 + D * log(2) + 
           log(TruncatedNormal::pmvnorm(mu_beta, Q_beta_inv, 
                                        lb = rep(0, D), ub = rep(Inf, D))[1])
    logml_truth[i] = LIL
    
    ep = list(H_k = Q_beta, H_k_inv = Q_beta_inv, m_k = mu_beta, 
              lambda = lambda)
    
    for (p in 1:N_REP) {
        
        samples = rtmvnorm(J, c(mu_beta), Q_beta_inv, rep(0, D), rep(Inf, D))
        u_df = preprocess(data.frame(samples), D, prior)
        
        approx[p, i] = hyb(u_df, ep)
        
        colnames(samples) = names(u_df)[1:D]
        lb = rep(0, D)
        ub = rep(Inf, D)
        names(lb) <- names(ub) <- colnames(samples)
        
        bridge_result <- bridge_sampler(samples = samples, 
                                        log_posterior = log_density,
                                        data = prior, lb = lb, ub = ub, 
                                        silent = TRUE)
        bridge[p, i] = bridge_result$logml
        
    }
    
    cat("D = ", D, ", ", "J = ", J, ";", "  ",
        "LIL = ",     format(round(logml_truth[i], 4),   nsmall = 4),
        " (hybrid: ", format(round(mean(approx[,i]), 4), nsmall = 4), 
        ", bridge: ", format(round(mean(bridge[,i]), 4), nsmall = 4), 
        ")", '\n', sep = "")
    
    i = i + 1
}


delta_hybrid = sweep(approx, 2, logml_truth) %>% colMeans
delta_bridge = sweep(bridge, 2, logml_truth) %>% colMeans

df = data.frame(hybrid = delta_hybrid, bridge = delta_bridge, 
                d = rep(D_vec, 2))

df_long = melt(df, id.vars = 'd', value.name = 'error', 
               variable.name = 'approx')

df_long %>% head

x11()
ggplot(df_long, aes(x = d, y = error, col = approx)) + geom_point(size = 2) + 
    geom_hline(yintercept = 0, col = 'red', size = 1) + 
    labs(y = "average error", title = "tMVN (Kappa = 4)")
    


df = data.frame(d = D_vec, 
                mae = colMeans(abs(delta)), 
                ae = colMeans(delta))

ggplot(df, aes(d, ae)) + geom_point(size = 2) + 
    geom_hline(yintercept = 0, col = 'red', size = 1) + 
    labs(y = "average error", title = "tMVN (Kappa = 4)")

















