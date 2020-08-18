


lik_vec = function(u, data) {
    beta = unname(unlist(u[1:p]))
    sigma2 = unname(unlist(u[D]))
    prod(dnorm(data$y, mean = data$X %*% beta, sd = sqrt(sigma2)))
}


data = list(y = c(y), X = X)

# samples from posterior stored row-wiseu
sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
beta_post = matrix(0, J, p)
for (j in 1:J) {
    beta_post[j,] = rmvnorm(1, mean = mu_star, sigma = sigmasq_post[j] * V_star)
}
u_samp = data.frame(beta_post, sigmasq_post)

# proposal distribution parameters
u_bar = unlist(unname(colMeans(u_samp))) # posterior mean
u_cov = unname(cov(u_samp))              # posterior covariance

# bridge sampling parameters
N_1 = J 
N_2 = J * 3
s_1 = N_1 / (N_1 + N_2)
s_2 = N_2 / (N_1 + N_2)

phat = 0 # starting point
phat_new = 10
tol = 1e-10
iter = 0
CONVERGE = FALSE

# set.seed(1)
while(!CONVERGE) {
    
    # form the numerator
    u_prop = rmvnorm(N_2, mean = u_bar, sigma = u_cov)    # sample from proposal
    g_prop = dmvnorm(u_prop, mean = u_bar, sigma = u_cov) # evaluate w/ proposal
    
    # evaluate samples from proposal using prior density
    ptheta_prop = dmvnorm(u_prop[,-D], mean = mu_beta, sigma = sigmasq * V_beta) * 
        MCMCpack::dinvgamma(c(u_prop[,D]), shape = a_0, scale = b_0)
    # evalute likelihood for each of the proposal samples
    lik_prop = apply(u_prop, 1, lik_vec, data = data)
    
    numer = lik_prop * ptheta_prop / (s_1 * lik_prop * ptheta_prop + s_2 * phat * g_prop)
    
    # form the denominator 
    g_post = dmvnorm(u_samp, mean = u_bar, sigma = u_cov)
    lik_post = apply(u_samp, 1, lik_vec, data = data)
    ptheta_post = dmvnorm(u_samp[,-D], mean = mu_beta, sigma = sigmasq * V_beta) * 
        MCMCpack::dinvgamma(c(u_samp[,D]), shape = a_0, scale = b_0)
    
    denom = g_post / (s_1 * lik_post * ptheta_post + s_2 * phat * g_post)
    
    phat_new = mean(numer) / mean(denom)
    
    iter = iter + 1    
    print(paste("iter =", iter))
    
    
    if (abs(log(phat) - log(phat_new)) < tol) {
        CONVERGE = TRUE
    }
    
    phat = phat_new
}

log(phat_new)
lil(y, X, prior, post)    # -256.7659


# ------------------------------------------------------------------------------
library(BayesianTools)
likelihood <- function(x) sum(msm::dtnorm(x, log = TRUE, lower = -1, upper = 1))

prior = createUniformPrior(lower = rep(-Inf, -Inf, 0), upper = rep(Inf, Inf,15))
bayesianSetup <- createBayesianSetup(likelihood = likelihood, prior = prior)
out = runMCMC(bayesianSetup = bayesianSetup, settings = list(iterations = 5000))

# plot(out)

# theoretical value
theory = log(1/(2^2))

marginalLikelihood(out)$ln.ML - theory


prior_density = function(u) {
    beta = unlist(unname(u[1:p]))
    sigmasq = unlist(unname(u[D]))
    dnorm(beta, mean = mu_beta, sd = sqrt(sigmasq)) * 
        dinvgamma(sigmasq, shape = a_0, scale = b_0)
} 


prior_sampler = function(n=1) {
   sigmasq = MCMCpack::rinvgamma(n, shape = a_0, scale = b_0)
   beta = rmvnorm(n, mean = mu_beta, sigma = sigmasq_post[j] * V_beta)
   cbind(beta), sigmasq)
}

prior = createPrior(density = prior_density, sampler = prior_sampler,
                    lower = c(rep(-10, p), 0),
                    upper = c(rep(20, p), 20))

likelihood = function(u) {
    beta = unlist(unname(u[1:p]))
    sigmasq = unlist(unname(u[D]))
    sum(dnorm(y, X %*% beta, sqrt(sigmasq), log = T))
}

mvnig_setup = createBayesianSetup(likelihood = likelihood, prior = prior)
settings = list(iterations = 1000)
out <- runMCMC(bayesianSetup = mvnig_setup, settings = settings)















