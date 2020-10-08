

library(BayesianTools)



# ------------------------------------------------------------------------------
library(ggplot2)

N = 50
D = 2


mu = 30
sigma_sq = 4
m_0 = 0
w_0 = 0.05
r_0 = 3
s_0 = 3

set.seed(123)
# generate 50 samples from N(mu, sigma_sq)
y = rnorm(N, mu, sqrt(sigma_sq))

prior = list(m_0 = m_0, w_0 = w_0, r_0 = r_0, s_0 = s_0, y = y)

ybar = mean(y)

# compute posterior parameters
m_n = (N * ybar + w_0 * m_0) / (N + w_0)
w_n = w_0 + N
r_n = r_0 + N
s_n = s_0 + sum((y - ybar)^2) + (N * w_0 / (N + w_0)) * (ybar - m_0)^2

p_y = pi^(-N / 2) * (w_0 / w_n)^(1/2) * gamma(r_n / 2) / gamma(r_0 / 2) * 
    s_0^(r_0 / 2) / s_n^(r_n / 2)

(LIL = log(p_y)) # -113.143

J = 1000
B = 100 # number of batch estimators



library(bridgesampling)
log_density = function(u, data) {
    -psi(u, data)
}





# load hml stuff ---------------------------------------------------------------
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")     
source("C:/Users/ericc/mlike_approx/paper_simulations/table1/normig_helper.R")

# compute the estimate for the i-th batch
hme    = numeric(B) 
hyb    = numeric(B) 
came   = numeric(B)
bridge = numeric(B)
set.seed(1)
for (i in 1:B) {
    
    # (0) sample from sigma_sq | y
    sigma_sq_post = MCMCpack::rinvgamma(J, shape = r_n / 2, scale = s_n / 2)
    # (1) sample from mu | sigma_sq, y
    mu_post = rnorm(J, m_n, sqrt(sigma_sq_post / w_n))
    
    # (1) compute hme ----------------------------------------------------------
    lik_j = numeric(J)
    for (j in 1:J) {
        # (2.1) compute the likelihood under the j-th posterior parameters
        lik_j[j] = 1 /  prod(dnorm(y, mu_post[j], sqrt(sigma_sq_post[j])))
    }
    # form the b-th HM estimator for the log integrated likelihood
    hme[i] = log(1 / (1 / J * sum(lik_j)))
    
    # (2) compute came ---------------------------------------------------------
    
    A_mu = c(min(mu_post), max(mu_post))
    A_sigmasq = c(min(sigma_sq_post), max(sigma_sq_post))
    
    # draw from the importance density N-IG
    sigmasq_mean = mean(sigma_sq_post)
    sigmasq_var = var(sigma_sq_post)
    
    # compute shape/scale parameters from posterior mean/var of sigmasq
    r_imp = sigmasq_mean^2/sigmasq_var + 2
    s_imp = sigmasq_mean * (r_imp - 1)
    
    sigmasq_s = MCMCpack::rinvgamma(J, shape = r_imp, scale = s_imp)
    mu_s      = rnorm(J, mean(mu_post), sqrt(mean(sigma_sq_post) / w_n))
    # mu_s      = rnorm(J, mean(mu_post), sqrt(mean(sigmasq_s) / w_n))
    
    # compute 1/s(theta) -- (K x 1) vector of evaluated densities
    # s_theta = dnorm(mu_s, m_n, sqrt(sigma_sq / w_n)) * 
    #     MCMCpack::dinvgamma(sigmasq_s, shape = r_n / 2, scale = s_n / 2)
    s_theta = dnorm(mu_s, mean(mu_post), sqrt(mean(sigma_sq_post) / w_n)) * 
        MCMCpack::dinvgamma(sigmasq_s, shape = r_imp, scale = s_imp)
    
    # compute prior density
    p_theta = dnorm(mu_s, m_0, sqrt(sigma_sq / w_0)) * 
        MCMCpack::dinvgamma(sigmasq_s, shape = r_0 / 2, scale = s_0 / 2)
    
    # compute likelihood
    lik_q = numeric(J)
    for (q in 1:J) {
        lik_q[q] = prod(dnorm(y, mean = mu_s[q], sd = sqrt(sigmasq_s[q])))
    }
    
    ind_A = (mu_s >=  A_mu[1] & mu_s <= A_mu[2]) & 
        (sigmasq_s >= A_sigmasq[1] & sigmasq_s <= A_sigmasq[2]) # 1_A (theta)
    
    came_approx = (1 / s_theta * lik_q * p_theta)[ind_A]
    came[i] = log(mean(came_approx))
    
    
    # (3) compute hybrid estimator ---------------------------------------------
    u_samps = data.frame(mu_post, sigma_sq_post)
    u_df = preprocess(u_samps, D, prior)
    
    # hybrid = hybrid_ml(D, u_df, J, prior)
    hybrid_v0 = hml_simple(D, u_df, J, prior)
    # hybrid_v0$param_out
    # hybrid_v0$zhat
    hyb[i] = hybrid_v0$zhat
    
    
    # (4) compute bridge estimator ---------------------------------------------
    lb <- c(-Inf, 0)
    ub <- c(Inf, Inf)
    u_samp = as.matrix(u_samps)
    colnames(u_samp) = names(u_df)[1:D]
    names(lb) <- names(ub) <- colnames(u_samp)
    bridge_result = bridge_sampler(samples = u_samp,
                                   log_posterior = log_density,
                                   data = prior, lb = lb, ub = ub, silent = TRUE)
    bridge[i] = bridge_result$logml
    
} # end main simulation loop

approx = data.frame(LIL, 
                    hme = hme[hme!=0],
                    hyb = hyb[hyb!=0],
                    bridge = bridge[bridge!=0],
                    came = came[came!=0])

error = data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
                   mae = colMeans(abs(LIL - approx)),
                   rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)
error




## arithmetic mean estimator --- pajor paper has this with 0 error ???
ame = numeric(B)
J = 10000
for (b in 1:B) {
    # (4) compute arithmetic mean estimator
    # sample from prior
    sigmasq_prior = MCMCpack::rinvgamma(J, shape = r_0 / 2, scale = s_0 / 2)
    mu_prior = rnorm(J, m_0, sqrt(sigmasq_prior / w_0))
    lik_prior = numeric(J)
    for (j in 1:J) {
        # (2.1) compute the likelihood under the j-th posterior parameters
        lik_prior[j] = prod(dnorm(y, mu_prior[j], sqrt(sigmasq_prior[j])))
    }
    ame[b] = log(mean(lik_prior))
}

## corrected arithmetic mean estimator
set.seed(123)
came = numeric(B)
J = 1000
for (b in 1:B) {
    
    # (1) sample from sigma_sq | y
    sigma_sq_post = MCMCpack::rinvgamma(J, shape = r_n / 2, scale = s_n / 2)
    # (0) sample from mu | sigma_sq, y
    mu_post = rnorm(J, m_n, sqrt(sigma_sq_post / w_n)) # (D x 1)
    
    A_mu = c(min(mu_post), max(mu_post))
    A_sigmasq = c(min(sigma_sq_post), max(sigma_sq_post))
    
    # draw from the importance density N-IG
    sigmasq_mean = mean(sigma_sq_post)
    sigmasq_var = var(sigma_sq_post)
    
    # compute shape/scale parameters from posterior mean/var of sigmasq
    r_imp = sigmasq_mean^2/sigmasq_var + 2
    s_imp = sigmasq_mean * (r_imp - 1)
    
    sigmasq_s = MCMCpack::rinvgamma(J, shape = r_imp, scale = s_imp)
    # mu_s      = rnorm(J, mean(mu_post), sqrt(mean(sigma_sq_post) / w_n))
    # mu_s      = rnorm(J, mean(mu_post), sqrt(mean(sigmasq_s) / w_n))
    
    # compute 1/s(theta) -- (K x 1) vector of evaluated densities
    s_theta = dnorm(mu_s, m_n, sqrt(sigma_sq / w_n)) * 
        MCMCpack::dinvgamma(sigmasq_s, shape = r_n / 2, scale = s_n / 2)
    
    # compute prior density
    p_theta = dnorm(mu_s, m_0, sqrt(sigma_sq / w_0)) * 
        MCMCpack::dinvgamma(sigmasq_s, shape = r_0 / 2, scale = s_0 / 2)
    
    # compute likelihood
    lik_q = numeric(J)
    for (q in 1:J) {
        lik_q[q] = prod(dnorm(y, mean = mu_s[q], sd = sqrt(sigmasq_s[q])))
    }
    
    ind_A = (mu_s >=  A_mu[1] & mu_s <= A_mu[2]) & 
        (sigmasq_s >= A_sigmasq[1] & sigmasq_s <= A_sigmasq[2]) # 1_A (theta)
    
    came_approx = (1 / s_theta * lik_q * p_theta)[ind_A]
    came[b] = log(mean(came_approx))
}


came %>% mean

approx_df = data.frame(mcmc = 1:B, hme = log(hme), hyb = hyb, came = came, ame)
approx_long = melt(approx_df, id.vars = 'mcmc')


ggplot(approx_long, aes(x = mcmc, y = value, col = variable)) + geom_point() +
    geom_hline(aes(yintercept = LIL), linetype = 'dashed', size = 0.9)


## compute average error, root mean squared error
mean(hyb)
mean(hme)
mean(ame)
mean(came)

mean(LIL - hyb)
mean(LIL - hme)
mean(LIL - ame)
mean(LIL - came)

sqrt(mean((LIL - hyb)^2))
sqrt(mean((LIL - hme)^2))
sqrt(mean((LIL - ame)^2))
sqrt(mean((LIL - came)^2))














# chib's implementation --------------------------------------------------------

#### using bayesiantools package

# prior density
density = function(par) {
    # normal-inverse gamma density
    d1 = dnorm(par[1], mean = m_n, sd = sqrt(par[2] / w_n))
    d2 = MCMCpack::dinvgamma(par[2], shape = r_n / 2, scale = s_n / 2)
    
    return(d1 * d2)
}


sampler = function(n = 1) {
    u2 = MCMCpack::rinvgamma(n, shape = r_0 / 2, scale = s_0 / 2)
    u1 = rnorm(n, m_0, sqrt(u2 / w_0))
    return(cbind(u1,u2))
}

prior <- createPrior(density = density, sampler = sampler, 
                     lower = c(-100,0), upper = c(100,500), best = NULL)


likelihood1 <- function(param){
    # pred = param[1] + param[2]*x + param[3] * x^2
    # singlelikelihoods = dnorm(y, mean = pred, sd = 1/(param[4]^2), log = TRUE)
    ll = dnorm(y, mean = param[1], sd = sqrt(param[2]), log = TRUE)
    return(sum(ll))  
}

setUp1 <- createBayesianSetup(likelihood1, prior = prior)
out1 <- runMCMC(bayesianSetup = setUp1)

out1$Z

M1 = marginalLikelihood(out1, numSamples = 5000, method = "Chib")
M1$ln.ML
LIL


u = getSample(out1)


# Use this prior in an MCMC 

ll <- function(x) sum(dnorm(x, log = TRUE)) # multivariate normal ll
bayesianSetup <- createBayesianSetup(likelihood = ll, prior = prior)

settings = list(iterations = 1000)
out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)















