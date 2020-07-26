


library(TruncatedNormal)
library(tmg)
library(mvtnorm)

# LEN_PATH  = "C:/Users/ericc/mlike_approx"
# setwd(LEN_PATH)
# 
# source("partition/partition.R")
# source("extractPartition.R")
# source("hybrid_approx.R")

setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")     
source("C:/Users/ericc/mlike_approx/truncate/regTN_helper.R")
D = 50
N = 200
I_D = diag(1, D)

n_samps = 10
J       = 1000
B       = 1000 # number of replications

hyb_fs  = numeric(B) # store harmonic mean estiator
hyb_ss  = numeric(B) # store harmonic mean estiator
hyb_ts  = numeric(B) # store harmonic mean estiator
hyb     = numeric(B) # store harmonic mean estiator


#### generate data -------------------------------------------------------------

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

(true_logml = lil_0 + D * log(2) + 
        log(TruncatedNormal::pmvnorm(mu_beta, Q_beta_inv, 
                                     lb = rep(0, D), ub = rep(Inf, D))[1]))



#### begin simulations ---------------------------------------------------------

for (b in 1:B) {
    
    # sample from posterior
    samples = data.frame(rtmvnorm(J, c(mu_beta), Q_beta_inv, rep(0, D), 
                                  rep(Inf, D)))

    u_df = preprocess(samples, D, prior)
    hml_approx = hml_const(1, D, u_df, J, prior)
    og_part = hml_approx$param_out %>% 
        dplyr::select(-c(psi_choice, logQ_cstar))
    ss_part = fit_resid(og_part, D, n_samps, prior)
    ts_part = fit_resid(ss_part, D, n_samps / 2, prior)
    
    hyb_fs[b] = hml_approx$const_vec
    hyb_ss[b] = log_sum_exp(unlist(compute_expterms(ss_part, D)))
    hyb_ts[b] = log_sum_exp(unlist(compute_expterms(ts_part, D)))
    
    hyb[b] = mean(c(hyb_fs[b], hyb_ss[b], hyb_ts[b]))
    
    #### (2) harmonic mean estimator
    # hme[b] = hme_approx(u_df, prior, J, D, N)
    
    #### (3) corrected arithmetic mean estimator (IS)
    # came[b] = came_approx(u_df, hml_approx, prior, post, J, D)
    
    if (b %% 100 == 0) { print(paste("iter:", b)) }
}





