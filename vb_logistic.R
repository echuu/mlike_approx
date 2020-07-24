
# number of MCMC samples
J = 4000
M = 9   # scale for posterior covariance



setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")     
source("C:/Users/ericc/mlike_approx/logistic/regTN_helper.R")

# ------------------------------------------------------------------------------

set.seed(1)
N = 25
D = 2
X = matrix(rnorm(N * D, 0, sqrt(1/N)), nrow = N)
beta = runif(D, 0.5, 3) # * sign(runif(D, -1,1)))
z = X %*% beta
pr = exp(z) / (1 + exp(z))
y = rbinom(N, 1, pr)

df = data.frame(y, X)

## prior params
tau = 1
prior = list(y = y, X = X, D = D, N = N, tau = tau)



# run VB to get posterior mean, covariance estimate
fit_vb <- vblogit(y, X, verb=TRUE) 
mu_beta = fit_vb$m # VB posterior mean
Sig_beta = matrix(fit_vb$S, D, D) # VB posterior covariance

## sample from the approximate posterior
u_samps = rmvnorm(J, mean = mu_beta, sigma = Sig_beta) %>% data.frame
# u_samps %>% head

# compute psi(u) for each of these
u_df = preprocess(u_samps, D, prior)
hml_approx = hml_const(1, D, u_df, J, prior)
# hml_approx$const_vec # -361.9583, -302.4611, -326.6722



## MCMC samples
bglm_fit = stan_glm(y ~ . -1, data = df, 
                    prior = normal(location = 0, scale = 1 / sqrt(tau)),
                    family = binomial(link = 'logit'),
                    refresh = 0)
summary(bglm_fit)

u_samps = bglm_fit %>% as.data.frame() # extract posterior samples

psi_u = apply(u_samps, 1, psi, prior = prior)
u_df = preprocess(u_samps, D, prior)

hml_approx = hml_const(1, D, u_df, J, prior) 


hml_approx$const_vec # -361.9583, -302.4611, -326.6722
hml_approx_vb$const_vec



G = 100

approx_vb     = numeric(G)
approx_vb_ss  = numeric(G)
approx_vb_ts  = numeric(G)

approx_mc     = numeric(G)
approx_mc_ss  = numeric(G)
approx_mc_ts  = numeric(G)

# true_logml    = numeric(G)

# d_0 = p / 2

# options(warn=0)
for (g in 1:G) {
    
    
    X = matrix(rnorm(N * D, 0, sqrt(1/N)), nrow = N)
    beta = runif(D, 0.5, 3) # * sign(runif(D, -1,1)))
    z = X %*% beta
    pr = exp(z) / (1 + exp(z))
    y = rbinom(N, 1, pr)
    
    ## prior params
    tau = 1
    prior = list(y = y, X = X, D = D, N = N, tau = tau)
    
    
    # run VB to get posterior mean, covariance estimate
    fit_vb <- vblogit(y, X, verb=TRUE) 
    mu_beta = fit_vb$m # VB posterior mean
    Sig_beta = matrix(fit_vb$S, D, D) # VB posterior covariance
    
    ## sample from the approximate posterior
    u_samps = rmvnorm(J, mean = mu_beta, sigma = Sig_beta) %>% data.frame

    # compute psi(u) for each of these
    u_df_vb = preprocess(u_samps, D, prior)
    hml_approx_vb = hml_const(1, D, u_df_vb, J, prior)
    # hml_approx$const_vec # -361.9583, -302.4611, -326.6722
    
    
    
    ## MCMC samples
    bglm_fit = stan_glm(y ~ . -1, data = df, 
                        prior = normal(location = 0, scale = 1 / sqrt(tau)),
                        family = binomial(link = 'logit'),
                        refresh = 0)
    u_samps = bglm_fit %>% as.data.frame() # extract posterior samples
    
    # psi_u = apply(u_samps, 1, psi, prior = prior)
    u_df_mcmc = preprocess(u_samps, D, prior)
    
    hml_approx_mcmc = hml_const(1, D, u_df_mcmc, J, prior) 
    
    approx_mc[g] = hml_approx_mcmc$const_vec
    approx_vb[g] = hml_approx_vb$const_vec  
    
    og_part = hml_approx_vb$param_out %>% 
        dplyr::select(-c(psi_choice, logQ_cstar))
    ss_fit = fit_resid(og_part, D, 5, prior)
    ts_fit = fit_resid(ss_fit, D, 5, prior)
    
    approx_vb_ss[g] = log_sum_exp(unlist(compute_expterms(ss_fit, D)))
    approx_vb_ts[g] = log_sum_exp(unlist(compute_expterms(ts_fit, D)))
    
    og_part = hml_approx_mcmc$param_out %>% 
        dplyr::select(-c(psi_choice, logQ_cstar))
    ss_fit = fit_resid(og_part, D, 5, prior)
    ts_fit = fit_resid(ss_fit, D, 5, prior)
    
    approx_mc_ss[g] = log_sum_exp(unlist(compute_expterms(ss_fit, D)))
    approx_mc_ts[g] = log_sum_exp(unlist(compute_expterms(ts_fit, D)))
    
    
    if (g %% 10 == 0) {
        print(paste("iter: ", g, " -- ", 
                    mean(c(approx_vb[g], approx_vb_ss[g], approx_vb_ts[g])),
                    ' (', 
                    mean(c(approx_mc[g], approx_mc_ss[g], approx_mc_ts[g])), 
                    ')', sep = ''))
    }
}


approx_vb     = numeric(G)
approx_vb_ss  = numeric(G)
approx_vb_ts  = numeric(G)

approx_mc     = numeric(G)
approx_mc_ss  = numeric(G)
approx_mc_ts  = numeric(G)

mean(approx_vb)
mean(approx_vb_ss)
mean(approx_vb_ts)

mean(approx_mc)
mean(approx_mc_ss)
mean(approx_mc_ts)

mean((approx_vb + approx_vb_ss + approx_vb_ts) / 3)
mean((approx_mc + approx_mc_ss + approx_mc_ts) / 3)

mean(approx_vb_ts)
mean(approx_mc_ts)

plot(approx_vb_ts, approx_mc_ts)
abline(0,1)





















### numerical integration ------------------------------------------------------

psi = function(u, prior) {
    
    y   = prior$y
    X   = prior$X
    D   = prior$D
    tau = prior$tau
    
    Xbeta = X %*% u
    
    loglik = (y * Xbeta - log(1 + exp(Xbeta))) %>% sum
    logprior = D / 2 * log(tau) - D / 2 * log(2 * pi) - 
        tau / 2 * sum(u * u)
    
    out = - loglik - logprior
    
    return(out)
    
} # end psi() function






fun <- function(x, y) exp(-n*x^2*y^4)
fun_logit = function(b1, b2) {
    
    # Xbeta = X %*% c(b1, b2)
    Xbeta = c(rep(b1, N) * X[,1] + rep(b2, N) * X[,2])
    Xbeta_0 = c(X %*% as.matrix(mu_beta))
    
    loglik = (y * (Xbeta - Xbeta_0) - log1pexp(Xbeta) +
                  log1pexp(Xbeta_0)) %>% sum
    
    logprior = D / 2 * log(tau) - D / 2 * log(2 * pi) - 
        tau / 2 * (b1^2 + b2^2)
    
    out = loglik + logprior
    exp(out)
}



Xbeta_0 = c(X %*% as.matrix(mu_beta))
loglik_0 = (y * Xbeta_0 - log(1 + exp(Xbeta_0))) %>% sum

exp(loglik_0) * fun_logit(u_df[1,1], u_df[1,2])
exp(-psi(unlist(unname(u_df[1,1:2])), prior))

library(pracma)
result = integral2(fun_logit, -2, 2, -5, 15, reltol = 1e-5)
result$Q
loglik_0 + log(result$Q) # -1.223014 for n = 1000
























### ----------------------------------------------------------------------------


# compare to glm() output
# fit_glm <- glm(y ~ -1+X, family=binomial)
# coefs <- cbind(vb=fit_vb$coef, glm=fit_glm$coef)
# plot(coefs, main="Estimates")
# abline(0,1)


## sample from the approximate posterior
u_samps = rmvnorm(J, mean = mu_beta, sigma = Sig_beta) %>% data.frame
# u_samps %>% head

# compute psi(u) for each of these
u_df = preprocess(u_samps, D, prior)

hml_approx = hml_const(1, D, u_df, J, prior)

hml_approx$const_vec # -361.9583, -302.4611, -326.6722
og_part = hml_approx$param_out %>% 
    dplyr::select(-c(psi_choice, logQ_cstar))

ss_part = fit_resid(og_part, D, 10, prior)
ts_part = fit_resid(ss_part, D, 10 / 2, prior)

log_sum_exp(unlist(compute_expterms(ss_part, D))) # -202.0286, -187.8912, -170.2579
log_sum_exp(unlist(compute_expterms(ts_part, D))) # -188.1801, -156.2409, -175.6461

mean(c(hml_approx$const_vec, 
       log_sum_exp(unlist(compute_expterms(ss_part, D))),
       log_sum_exp(unlist(compute_expterms(ts_part, D))))) # -250.7224, -215.531, -233.7635

length(unlist(compute_expterms(ts_part, D)))

# ------------------------------------------------------------------------------
df = data.frame(y, X)

bglm_fit = stan_glm(y ~ . -1, data = df, 
                    prior = normal(location = 0, scale = 1 / sqrt(tau)),
                    family = binomial(link = 'logit'),
                    refresh = 0)

# summary(bglm_fit)

u_samps_post = bglm_fit %>% as.data.frame() # extract posterior samples
u_df_post = preprocess(u_samps_post, D, prior)
hml_approx_post = hml_const(1, D, u_df_post, J, NULL)
hml_approx_post$const_vec # -379.2573, -311.5547, -393.1815

og_part_post = hml_approx_post$param_out %>% 
    dplyr::select(-c(psi_choice, logQ_cstar))

ss_part_post = fit_resid(og_part_post, D, 10, prior)
ts_part_post = fit_resid(ss_part_post, D, 10 / 2, prior)

length(unlist(compute_expterms(ts_part_post, D)))

log_sum_exp(unlist(compute_expterms(ss_part_post, D))) # -197.8441, -176.446, -201.5819
log_sum_exp(unlist(compute_expterms(ts_part_post, D))) # -164.5808, -170.1808, -180.243

mean(c(hml_approx_post$const_vec, 
       log_sum_exp(unlist(compute_expterms(ss_part_post, D))),
       log_sum_exp(unlist(compute_expterms(ts_part_post, D))))) # -247.2274, -219.3938, -238.358






