
# hme_estimator.R


# -----------------------------------------------------------------------------
library('dplyr')
library('ggplot2')
library('MCMCpack') # for rinvgamma() function
library('tree')     # plotting partitions of a fitted decision tree
library('reshape2')

# lenk paper simulation 

mu = 30
sigma_sq = 4
m_0 = 0
w_0 = 0.05
r_0 = 3
s_0 = 3

N = 100
set.seed(123)
# generate 50 samples from N(mu, sigma_sq)
y = rnorm(N, mu, sqrt(sigma_sq))

ybar = mean(y)

# compute posterior parameters
m_n = (N * ybar + w_0 * m_0) / (N + w_0)
w_n = w_0 + N
r_n = r_0 + N
s_n = s_0 + sum((y - ybar)^2) + (N * w_0 / (N + w_0)) * (ybar - m_0)^2

p_y = pi^(-N / 2) * (w_0 / w_n)^(1/2) * gamma(r_n / 2) / gamma(r_0 / 2) * 
    s_0^(r_0 / 2) / s_n^(r_n / 2)

LIL = log(p_y) # -214.5399 (paper says -117, but difference arises from RNG)


# generate samples from the posterior probability to form the HME estimator

J = 40 # number of random draws used per estimate
B = 20 # number of batch estimators


# compute the estimate for the i-th batch
lil_hat = numeric(B) # store the log integrated likelihood for each batch
for (b in 1:B) {
    
    # (0) sample from mu | sigma_sq, y
    mu_post = rnorm(J, m_n, sqrt(sigma_sq / w_n)) # (D x 1)
    
    # (1) sample from sigma_sq | y
    sigma_sq_post = MCMCpack::rinvgamma(J, shape = r_n / 2, scale = s_n / 2)
    
    # (2) compute the harmonic mean estimator
    lik_j = numeric(J)
    for (j in 1:J) {
        # (2.1) compute the likelihood under the j-th posterior parameters
        lik_j[j] = 1 /  prod(dnorm(y, mu_post[j], sqrt(sigma_sq_post[j])))
    }
    
    # form the b-th HM estimator for the log integrated likelihood
    lil_hat[b] = 1 / (1 / J * sum(lik_j))
    
}

hme_df = data.frame(mcmc = 1:B, hme = log(lil_hat), lil = LIL)

# generate figure 1 in lenk 2009
ggplot(hme_df, aes(x = mcmc, y = hme)) + geom_point() +
    geom_hline(aes(yintercept = LIL), linetype = 'dashed', size = 0.9)


# corrected arithmetic mean estimator ------------------------------------------



log_s_theta = function(mu, sigmasq, m, w, a, b) {
    
    out = a * log(b) - d / 2 * log(2 * pi) - 0.5 * log_det(V_beta) - lgamma(a) -
        (a + d / 2 + 1) * log(sigmasq) - 
        1 / sigmasq * (b + 0.5 * t(beta - mu_beta) %*% V_beta_inv %*% 
                           (beta - mu_beta))
    
    
}




# store log marginal likelihood estimators for B batches

lil_came = numeric(B) # CAME estimator
lil_hme  = numeric(B) # store the log integrated likelihood for each batch

J = 40 # number of draws from the importance function
B = 20 # number of batch estimators

for (b in 1:B) {
    
    # (0) sample from mu | sigma_sq, y
    mu_post = rnorm(J, m_n, sqrt(sigma_sq / w_n)) # (D x 1)
    
    # (1) sample from sigma_sq | y
    sigma_sq_post = MCMCpack::rinvgamma(J, shape = r_n / 2, scale = s_n / 2)
    
    # (2) compute corrected arithmetic mean estimator  -------------------------
    A_mu = c(min(mu_post), max(mu_post))
    A_sigmasq = c(min(sigma_sq_post), max(sigma_sq_post))
    
    # draw from the importance density N-IG
    mu_s      = rnorm(K, m_n, sqrt(sigma_sq / w_n))
    sigmasq_s = MCMCpack::rinvgamma(K, shape = r_n / 2, scale = s_n / 2)
    
    # compute 1/s(theta) -- (K x 1) vector of evaluated densities
    s_theta = dnorm(mu_s, m_n, sqrt(sigma_sq / w_n)) * 
        MCMCpack::dinvgamma(sigmasq_s, shape = r_n / 2, scale = s_n / 2)
    
    s_tilde = log_mvnig(c(mu_s, s))
    
    
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
    
    lil_hat = 1 / s_theta * lik_q * p_theta * ind_A
    # --------------------------------------------------------------------------
    
    # (2) compute the harmonic mean estimator
    lik_j = numeric(J)
    for (j in 1:J) {
        # (2.1) compute the likelihood under the j-th posterior parameters
        lik_j[j] = 1 /  prod(dnorm(y, mu_post[j], sqrt(sigma_sq_post[j])))
    }
    
    # form the b-th HM estimator for the log integrated likelihood
    lil_hme[b] = log(1 / (1 / J * sum(lik_j)))
    lil_came[b] = log(mean(lil_hat))
    
}

hme_df = data.frame(mcmc = 1:B, hme = lil_hme, came = lil_came, lil = LIL)

hme_df_long = melt(hme_df, id.vars = "mcmc")

ggplot(hme_df_long, aes(x = mcmc, y = value, col = variable)) + geom_point() +
    geom_hline(aes(yintercept = LIL), linetype = 'dashed', size = 0.9)































