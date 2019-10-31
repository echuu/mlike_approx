
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

LIL = log(p_y) # -118.332 (paper says -117, but difference arises from RNG)

# define the function psi()
psi = function(y, mu, sigma_sq, m_0, w_0, r_0, s_0) {
    like = dnorm(y, mu, sqrt(sigma_sq))
    p_mu = dnorm(mu, m_0, sqrt(sigma_sq / w_0))
    p_sigma_sq = dinvgamma(sigma_sq, r_0 / 2, s_0 / 2)
    return(log(like * p_mu * p_sigma_sq))
}

# generate samples from the posterior probability to form the HME estimator

J = 1000 # number of random draws used per estimate

# (0) sample from mu | sigma_sq, y
mu_post = rnorm(J, m_n, sqrt(sigma_sq / w_n)) # (D x 1)

# (1) sample from sigma_sq | y
sigma_sq_post = rinvgamma(J, shape = r_n / 2, scale = s_n / 2)

# (2) for each u drawn from the posterior, evaluate psi(u) = psi(mu, sigma_sq)
psi_u = psi(y, mu_post, sigma_sq_post, m_0, w_0, r_0, s_0)
u_df = data.frame(mu = mu_post, sigsq = sigma_sq_post, psi_u = psi_u)

# (3) fit decision tree
u_tree = tree(psi_u ~ mu + sigsq, u_df)

# (4) plot the partition over the parameter space
partition.tree(u_tree, cex = 1)

# (5) look at posterior distributions of mu, sigma_sq to see if there are more
#     partitions in areas of higher posterior probability

hist(mu_post)
hist(sigma_sq_post)


par(mfrow = c(3, 1))





