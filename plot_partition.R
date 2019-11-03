

library('dplyr')
library('ggplot2')
library('MCMCpack')  # for rinvgamma() function
library('tree')      # plotting partitions of a fitted decision tree
library('mvtnorm')   # multivariate normal density

set.seed(123)

N = 50

# y | mu ~ N (mu, Sigma)
mu = c(3, 1)
d = length(mu)
Omega = matrix(c(1, 0.7, 0.7, 1), nrow = 2)
xi_0 = d / N
Sigma = xi_0 * Omega

# mu ~ N (mu_0, Sigma_0)
mu_0 = c(0, 0)
Sigma_0 = diag(1, 2)

# generate N samples from y | mu
y = mvrnorm(N, mu, Sigma)
ybar = apply(y, 2, mean) # mean of each component

# define the posterior distribution parameters: mu | y ~ N (mu_n, Sigma_n)
Sigma_N = solve(solve(Sigma_0) + N * solve(Sigma))
mu_N = Sigma_N %*% (solve(Sigma_0) %*% mu_0 + N * solve(Sigma) %*% ybar)

# generate samples from the posterior probability to form the HME estimator
J = 1000 # number of random draws used per estimate

# (0) sample from mu | sigma_sq, y
mu_post = mvrnorm(J, mu_N, Sigma_N) # (D x 1)

# plot the points drawn from the posterior distribution
mu_df = data.frame(mu1 = mu_post[,1], mu2 = mu_post[,2])
ggplot(mu_df, aes(mu1, mu2)) + geom_point() + theme_bw()


# define function psi()
psi = function(y, mu, Sigma, mu_0, Sigma_0) {
    log_like = dmvnorm(y, mu, Sigma) # evaluate y row-wise
    
    # likelihood needs to be computed for each of the parameter settings 
    
    p_mu = dmvnorm(mu, mu_0, Sigma_0)
    
    return(list(psi_u = log(like * p_mu * p_sigma_sq),
                p_mu = p_mu, p_sigma_sq = p_sigma_sq))
}


# TODO -------------------------------------------------------------------------

# (3) fit decision tree
u_tree = tree(psi_u ~ mu + sigsq, u_df)

# (4) plot the partition over the parameter space
partition.tree(u_tree, cex = 1)

# (5) look at posterior distributions of mu, sigma_sq to see if there are more
#     partitions in areas of higher posterior probability





