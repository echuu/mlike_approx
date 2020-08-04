
library(bridgesampling)

# proposal = 'normal'


stancodeH1 <- 'data {
  int<lower=1> n;    // number of observations
  int<lower=1> p;    // covariate dimension
  vector[n] y;       // observations
  matrix[n, p] X;    // design matrix
  vector[p] mu_0;    // prior mean
  matrix[p, p] I_p;  // (pxp) identity
  matrix[n, n] I_n;  // (nxn) identity
  real<lower=0> a_0; // prior shape
  real<lower=0> b_0; // prior scale
}
parameters {
  real<lower=0> sigmasq; // variance
  vector[p] beta;        // coefficient vector
}
model {
  target += inv_gamma_lpdf(sigmasq | a_0, b_0);
  target += multi_normal_lpdf(beta | mu_0, sigmasq * I_p);
  target += multi_normal_lpdf(y | X * beta, sigmasq * I_n);
}
'

# compile models
stanfitH1 <- sampling(stanmodelH1, data = list(n = N, p = p,
                                               y = c(y), X = X, 
                                               mu_0 = mu_beta, 
                                               I_p = I_p, I_n = I_N,
                                               a_0 = a_0, b_0 = b_0),
                      iter = 2000, warmup = 1000, chains = 3, cores = 1)

post_params = extract(stanfitH1, pars = c("beta", "sigmasq"))
head(post_params$beta)

plot(post_params$beta)
plot(u_df[,-(3:4)])


# compute log marginal likelihood via bridge sampling for H1
H1.bridge = bridge_sampler(stanfitH1, silent = FALSE)
H1.bridge$logml
H1.bridge$niter

u_df %>% dim


# samples <- rmvnorm(1e4, mean = rep(0, 2), sigma = diag(2))
# colnames(samples) <- c("x1", "x2")


# this function is just -psi(u)
log_density = function(u, data) {
    
    beta = unname(unlist(u[1:p]))
    sigma2 = unname(unlist(u[D]))
    
    sum(dnorm(y, mean = X %*% beta, sd = sqrt(sigmasq), log = T)) + 
        c(a_0 * log(b_0) - p / 2 * log(2 * pi) - 
              0.5 * log_det(V_beta) - lgamma(a_0) -
              (a_0 + p / 2 + 1) * log(sigmasq) - 
              1 / sigmasq * (b_0 + 0.5 * t(beta - mu_beta) %*% V_beta_inv %*% 
                                 (beta - mu_beta)))
}

log_density = function(u, data) {
    -psi(u, data)
}


# sample from posterior
sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
beta_post = matrix(0, J, p)
for (j in 1:J) {
    beta_post[j,] = rmvnorm(1, mean = mu_star, sigma = sigmasq_post[j] * V_star)
}
u_samp = data.frame(beta_post, sigmasq_post)
u_df = preprocess(u_samp, D, prior)

u_samp = as.matrix(u_samp)
colnames(u_samp) = names(u_df)[1:D]

lb <- c(rep(-Inf, p), 0)
ub <- c(rep(Inf, p), Inf)
names(lb) <- names(ub) <- colnames(u_samp)

bridge_result <- bridge_sampler(samples = u_samp, log_posterior = log_density,
                                data = prior, lb = lb, ub = ub, silent = TRUE)

# bridge_result_warp <- bridge_sampler(samples = u_samp, 
#                                      log_posterior = log_density,
#                                      method = 'warp3', data = NULL, 
#                                      lb = lb, ub = ub, silent = TRUE)


# bridge_result$method
# bridge_result$niter
# bridge_result_warp$logml
bridge_result$logml
lil(y, X, prior, post)    # -256.7659

























