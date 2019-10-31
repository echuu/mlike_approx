
# hme_estimator.R


# -----------------------------------------------------------------------------
library('dplyr')
library('ggplot2')
library('MCMCpack') # for rinvgamma() function
library('tree')     # plotting partitions of a fitted decision tree

# lenk paper simulation 

mu = 30
sigma_sq = 4
m_0 = 0
w_0 = 0.05
r_0 = 3
s_0 = 3

N = 50

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

LIL = log(p_y) # -118.332 (paper says -117, but difference arises from RNG)


# generate samples from the posterior probability to form the HME estimator

J = 1000 # number of random draws used per estimate
B = 1000 # number of batch estimators


# compute the estimate for the i-th batch
lil_hat = numeric(B) # store the log integrated likelihood for each batch
for (b in 1:B) {
    
    # (0) sample from mu | sigma_sq, y
    mu_post = rnorm(J, m_n, sqrt(sigma_sq / w_n)) # (D x 1)
    
    # (1) sample from sigma_sq | y
    sigma_sq_post = rinvgamma(J, shape = r_n / 2, scale = s_n / 2)
    
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