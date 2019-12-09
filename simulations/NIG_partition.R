
# NIG_partition.R --------------------------------------------------------------

# normal-inverse gamma example

# model:  y_i | mu, sigma_sq ~ N   (mu, sigma_sq)
#         mu  | sigma_sq     ~ N   (m_0, sigma_sq / w_0)
#               sigma_sq     ~ IG  (r_0 / 2, s_0 / 2)

# posterior(s):  mu        | sigma_sq, y ~ N   (m_n, sigma_sq / w_n)
#                sigma_sq  | y           ~ IG  (r_n / 2, s_n / 2)


# additional notes: requires latest version of R to load the 'tree' package
# uncomment lines to install missing packages: 
#     install.packages("dplyr")
#     install.packages("ggplot2")
#     install.packages("MCMCpack")
#     install.packages("tree")

# -----------------------------------------------------------------------------
library('dplyr')
library('ggplot2')
library('MCMCpack')  # for rinvgamma() function
library('tree')      # plotting partitions of a fitted decision tree

library('invgamma')  # dinvgamma(x, shape, scale)


f = function(x, shape = 1/2, scale = 1/2) {
    scale^shape / gamma(shape) * x^(- shape - 1) * exp(-scale / x)
}



set.seed(123)

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

LIL = log(p_y) # -113.143 (paper says -117, but difference arises from RNG)

# y        : (N x 1)
# mu       : (J x 1) J samples of mu from the posterior
# sigma_sq : (J x 1) J samples of sigma^2 from the posterior
# m_0      : (1 x 1) prior mean
# w_0      : (1 x 1) scale the prior variance
# r_0      : (1 x 1) prior shape
# s_0      : (1 x 1) prior scale
psi_p = function(y, mu, sigma_sq, m_0, w_0, r_0, s_0) {
    
    # for  y_1:N, evaluate the likelihood for j-th set of parameter values
    loglik = numeric(length(mu))
    
    for (j in 1:length(mu)) {
        loglik[j] = sum(dnorm(y, mu[j], sqrt(sigma_sq[j]), log = TRUE))
    }
    
    # following two computations are unused for now
    log_p_mu = dnorm(mu, m_0, sqrt(sigma_sq / w_0), log = TRUE)
    log_p_sigma_sq = invgamma::dinvgamma(sigma_sq, shape = r_0 / 2, 
                                         scale = s_0 / 2, log = TRUE)
    
    return(loglik) 
} # end psi() function
# ------------------------------------------------------------------------------

# generate samples from the posterior probability to form the HME estimator
J = 1000 # number of random draws used per estimate

# (0) sample from mu | sigma_sq, y
mu_post = rnorm(J, m_n, sqrt(sigma_sq / w_n)) # (D x 1)

# (1) sample from sigma_sq | y
sigma_sq_post = MCMCpack::rinvgamma(J, shape = r_n / 2, scale = s_n / 2)

# (2) for each u drawn from the posterior, evaluate psi(u) = psi(mu, sigma_sq)
psi_u = psi_p(y, mu_post, sigma_sq_post, m_0, w_0, r_0, s_0) # (J x 1)
u_df = data.frame(mu = mu_post, sigsq = sigma_sq_post, psi_u = psi_u) # (J x 3)

# (3) fit decision tree
u_tree = tree(psi_u ~ mu + sigsq, u_df)

# (4) plot the partition over the parameter space
par(mfrow = c(1,1))
plot(u_df[,1], u_df[,2], pch = 20, cex = 0.9, col = "pink",
     xlab = 'mu', ylab = 'sigma_sq', main = 'N = 1000, J = 5000')
partition.tree(u_tree, add = TRUE, cex = 0.01, ordvars = c("mu", "sigsq"))

partition.tree(u_tree, add = F, cex = 0.01, ordvars = c("mu", "sigsq"))


# ------------------------------------------------------------------------------

# par(mfrow = c(1,2))
# hist(mu_post)
# hist(sigma_sq_post)
# par(mfrow = c(1,1))


# partition.tree(u_tree, add = F, cex = 0.01, ordvars = c("mu", "sigsq"))

# par(mfrow = c(1, 3))





