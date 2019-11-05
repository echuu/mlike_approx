

library('dplyr')
library('ggplot2')
library('MCMCpack')  # for rinvgamma() function
library('tree')      # plotting partitions of a fitted decision tree
library('mvtnorm')   # multivariate normal density

set.seed(123)

# N = 50, 100, 1000
N = 100

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

# define function psi()
psi = function(y, mu, Sigma, mu_0, Sigma_0) {
    
    # for the vectors y, evaluate the log-like for j-th set of parameter values
    loglik = numeric(nrow(mu))
    
    # likelihood needs to be computed for each of the parameter settings 
    for (j in 1:nrow(mu)) {
        loglik[j] = sum(dmvnorm(y, mu[j,], Sigma, log = TRUE))
    }
    
    # log prior of mu -- unused for now
    log_p_mu = dmvnorm(mu, mu_0, Sigma_0, log = TRUE)
    
    return(loglik)
} # end psi() function

# generate samples from the posterior probability to form the HME estimator
# J = 500, 1000, 5000
J = 5000 # number of random draws used per estimate

# (0) sample from mu | y
mu_post = mvrnorm(J, mu_N, Sigma_N) # (D x 1)

# plot the points drawn from the posterior distribution
mu_df = data.frame(mu1 = mu_post[,1], mu2 = mu_post[,2])
ggplot(mu_df, aes(mu2, mu1)) + geom_point() + theme_bw()

# (1) evaluate psi(u) = psi(mu)
psi_u = psi(y, mu_post, Sigma, mu_0, Sigma_0)
u_df = data.frame(mu1 = mu_post[,1], mu2 = mu_post[,2], psi_u = psi_u) # (J x 3)

# (2) fit decision tree
u_tree = tree(psi_u ~ mu1 + mu2, u_df)

# (3) plot the partition over the parameter space
partition.tree(u_tree, cex = 1, ordvars = c("mu1", "mu2"))


par(mfrow = c(1,2))
hist(mu_post[,1])
hist(mu_post[,2])


# overlay partition on scatterplot of the posterior distribution
par(mfrow = c(1,1))
plot(mu_df[,1], mu_df[,2], pch = 20, cex = 0.8, col = "cyan",
     xlab = 'mu1', ylab = 'mu2', main = 'N = 100, J = 5000')
partition.tree(u_tree, add = TRUE, cex = 0.01, ordvars = c("mu1", "mu2"))






