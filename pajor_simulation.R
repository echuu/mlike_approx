
# pajor_simulation.R -- implement the marginal likelihood estimation method
# described in arithmetic mean paper



# 3.1 Conjugate Normal Model
# y_t | mu, sigma^2 ~ N  (y_t     | mu, sigma^2), t = 1,...,N
#  mu | sigma^2     ~ N  (mu      | m_0, sigma^2 / w_0)
#       sigma^2     ~ IG (sigma^2 | r_0 / 2, s_0 / 2)

# conditional posterior distribution for mu is normal: 
# 



N = 100 # number of observations generated from the above model
MC_iter = 1e4 # number of monte carlo simualations 


# 10,000 iterations in each MC procedure
#     1,000 re-estimations for N = 100 observations were carried out
#     (normal with mean 0, variance 4)

# normal params for N observations
mu = 0
sigma2 = 4

# prior hyperparameters --------------------------------------------------------

# normal params for mu
m_0 = 0
w_0 = 0.05

# inv-gamma params for sigma^2
r_0 = 3
s_0 = 3


# generate y1,...,y_N using distribution specified above -----------------------
y = rnorm(N, mean = mu, sd = sqrt(sigma2))



# harmonic mean estimator


# arithmetic mean estimator