
# NIX_partition.R --------------------------------------------------------------

# install.packages("LaplacesDemon)
# install.packages("metRology)

library("LaplacesDemon")  # sample form inverse chi square distribution
library("metRology")      # sample from scaled, shifted t distribution
library("tree")

# normal-inverse-chisquare example

# y_i | mu, sigma_sq ~ N (mu, sigma_sq)
# mu  | sigma_sq     ~ N (mu_0, sigma_sq / k_0)
#       sigma_sq     ~ invchisq (nu_0, sigma_sq_0)

set.seed(123)

N = 100
mu = 30
sigma_sq = 4
mu_0 = 0
k_0 = 0.05
nu_0 = 5        # prior df for inverse chisquare
sigma_sq_0 = 1  # prior scale for inverse chisquare

# generate N samples from N(mu, sigma_sq)
y = rnorm(N, mu, sqrt(sigma_sq))  # (N x 1)

# intermediate quantities
ybar = mean(y)
s_y = sum((y - ybar)^2)

# define the psi() function used to evaluate the posterior samples -------------
psi = function(y, mu, sigma_sq) {
    
    # for  y_1:N, evaluate the likelihood for j-th set of parameter values
    loglik = numeric(length(mu))
    
    for (j in 1:length(mu)) {
        # evaluate the log-likelihood of the j-th posterior sample
        loglik[j] = sum(dnorm(y, mu[j], sqrt(sigma_sq[j]), log = TRUE))
    }
    
    return(loglik)
} # end psi() function
# ------------------------------------------------------------------------------

J = 5000

# define the posterior parameter estimates
nu_N = nu_0 + N
sigma_sq_N = 1 / nu_N * (nu_0 * sigma_sq_0 + s_y + 
                             N * k_0 / (k_0 + N) * (mu_0 - ybar)^2)
k_N = k_0 + N
mu_N = (k_0 * mu_0 + N * ybar) / k_N

## draw from posterior distributions

# sample from mu | y ~ t_{nu_N} (mu_N, sigma_sq_N / k_N)
mu_post = rt.scaled(J, df = nu_N, mu_N, sqrt(sigma_sq_N / k_N)) 
# hist(mu_post)

# posterior of sigma_sq | y ~ invchisq(nu_n, sigma_sq_n)
sigma_sq_post = rinvchisq(J, df = nu_N, scale = sigma_sq_N)
# hist(sigma_sq_post) # looks about right

# (1) evaluate psi(u) = psi(mu)
psi_u = psi(y, mu_post, sigma_sq_post)
u_df = data.frame(mu = mu_post, sigma_sq = sigma_sq_post, psi_u = psi_u)

# (2) fit decision tree
u_tree = tree(psi_u ~ mu + sigma_sq, u_df)

# (3) plot the partition over the parameter space
# partition.tree(u_tree, cex = 1)



# compare the partitioned parameter space to the posterior distributions
# par(mfrow = c(1,2))
# hist(mu_post)
# hist(sigma_sq_post)


# overlay partition on scatterplot of the posterior distribution
par(mfrow = c(1,1))
plot(u_df[,2], u_df[,1], pch = 20, cex = 0.9, col = "pink",
     xlab = 'sigma_sq', ylab = 'mu', main = 'N = 1000, J = 5000')
partition.tree(u_tree, add = TRUE, cex = 0.01)


# overlay partitions corresponding to different # of MC samples ----------------











