

library(BayesianTools)


# Creating test data with quadratic relationship
sampleSize = 30
x <- (-(sampleSize-1)/2):((sampleSize-1)/2)
y <-  1 * x + 1*x^2 + rnorm(n=sampleSize,mean=0,sd=10)
# plot(x,y, main="Test Data")

# likelihoods for linear and quadratic model 
likelihood1 <- function(param){
    pred = param[1] + param[2]*x + param[3] * x^2
    singlelikelihoods = dnorm(y, mean = pred, sd = 1/(param[4]^2), log = TRUE)
    return(sum(singlelikelihoods))  
}
likelihood2 <- function(param){
    pred = param[1] + param[2]*x 
    singlelikelihoods = dnorm(y, mean = pred, sd = 1/(param[3]^2), log = TRUE)
    return(sum(singlelikelihoods))  
}

setUp1 <- createBayesianSetup(likelihood1, 
                              lower = c(-5,-5,-5,0.01), 
                              upper = c(5,5,5,30))
setUp2 <- createBayesianSetup(likelihood2, 
                              lower = c(-5,-5,0.01), 
                              upper = c(5,5,30))

out1 <- runMCMC(bayesianSetup = setUp1)
M1 = marginalLikelihood(out1, start = 1000)

out2 <- runMCMC(bayesianSetup = setUp2)
M2 = marginalLikelihood(out2, start = 1000)

M1$ln.ML

M2$ln.ML

# ------------------------------------------------------------------------------
library(ggplot2)

N = 50
D = 2


mu = 30
sigma_sq = 4
m_0 = 0
w_0 = 0.05
r_0 = 3
s_0 = 3

prior = list(m_0 = m_0, w_0 = w_0, r_0 = r_0, s_0 = s_0, y = y)

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

LIL = log(p_y) # -113.143

J = 1000
B = 1000 # number of batch estimators


# load hml stuff ---------------------------------------------------------------
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")     
source("C:/Users/ericc/mlike_approx/paper_simulations/table1/normig_helper.R")


# (1) sample from sigma_sq | y
sigma_sq_post = MCMCpack::rinvgamma(J, shape = r_n / 2, scale = s_n / 2)
mu_post = rnorm(J, m_n, sqrt(sigma_sq_post / w_n)) # (D x 1)

u_samps = data.frame(mu_post, sigma_sq_post)
u_df = preprocess(u_samps, D, prior)

hml_approx = hml_const(1, D, u_df, J, prior)
hml_approx$const_vec
hml_approx$param_out


# compute the estimate for the i-th batch
hme = numeric(B) # store the log integrated likelihood for each batch
hyb = numeric(B) # store the log integrated likelihood for each batch

for (b in 1:B) {
    
    # (0) sample from sigma_sq | y
    sigma_sq_post = MCMCpack::rinvgamma(J, shape = r_n / 2, scale = s_n / 2)
    # (1) sample from mu | sigma_sq, y
    mu_post = rnorm(J, m_n, sqrt(sigma_sq_post / w_n))
    
    # (2) compute the harmonic mean estimator
    lik_j = numeric(J)
    for (j in 1:J) {
        # (2.1) compute the likelihood under the j-th posterior parameters
        lik_j[j] = 1 /  prod(dnorm(y, mu_post[j], sqrt(sigma_sq_post[j])))
    }
    # form the b-th HM estimator for the log integrated likelihood
    hme[b] = 1 / (1 / J * sum(lik_j))
    
    # (3) compute hybrid estimator
    u_samps = data.frame(mu_post, sigma_sq_post)
    u_df = preprocess(u_samps, D, prior)
    
    hml_approx = hml_const(1, D, u_df, J, prior)
    
    hyb[b] = hml_approx$const_vec
    
}

approx_df = data.frame(mcmc = 1:B, hme = log(hme), hyb = hyb, lil = LIL)
approx_long = melt(approx_df, id.vars = 'mcmc')


ggplot(approx_long, aes(x = mcmc, y = value, col = variable)) + geom_point() +
    geom_hline(aes(yintercept = LIL), linetype = 'dashed', size = 0.9)


## compute average error, root mean squared error
mean(hyb)
mean(log(hme))

mean(LIL - log(hme))
mean(LIL - hyb)

sqrt(mean((LIL - log(hme))^2))
sqrt(mean((LIL - hyb)^2))





















