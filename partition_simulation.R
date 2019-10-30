

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








