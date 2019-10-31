

## preliminary steps to check if the partitioning is doing the right thing

# (i)  simulate data drawn from the posterior distribution, evaluate it 
#      to get the (u_1, psi(u_1)),...,(u_T, psi(u_T)) pairs
# (ii) feed the data into decision tree to obtain the partition on the 
#      parameter space

# define the negative log posterior -- from the setup of the problem, we assume
# that we can evaluate this function

# psi() -- evaluate the d-dimensional input vector u
# inputs: 
#         u     : the input vector whose density is to be calculated
#         d     : dimension of input
#         mu    : posterior mean
#         Omega : scaled covariance matrix
#         N     : N > d, serves as the 'ficitious proxy for sample size'
psi = function(u, d, mu, Omega, N) {
    
    # scale the (fixed) covariance matrix to make the target distribution
    # relatively concentrated
    Sigma = d / N * Omega  # (d x d) covariance matrix
    
    return(-0.5 * t(u - mu) %*% solve(Sigma) %*% (u - mu))
} # end psi() function


psi = function(u, mu, sigma) {
    
    return(-0.5 * (u - mu)^2 / sigma^2)
} # end psi() function

# samples from gamma = N( u | mu, Sigma), d-dim gaussian

d = 1
N = 100
mu = rep(0, d)
sigma = d / N * 1


post_sample = rnorm(N, mu, sd = sqrt(d/N*1))
psi_sample = psi(post_sample, mu, sqrt(sigma))

df_post = data.frame(x = post_sample, y = psi_sample)

# extract the decision boundaries from the regression tree
reg.tree = rpart(y ~ x, data = df_post)
rpart.plot(reg.tree, type = 4)

reg.tree$splits

# -----------------------------------------------------------------------------












