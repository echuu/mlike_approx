


## simulation code 



D = c(4) # test for smalller dimensions for now
N = c(100) # for testing -- comment this line to perform ext. analysis

## algorithm settings ----------------------------------------------------------
J         = 2000         # number of MC samples per approximation
# ------------------------------------------------------------------------------

set.seed(1)
D_MAX = 100
D_vec = seq(10, D_MAX, 2)
N_REP = 50
KAPPA = 4 # (# MCMC samples) / D = J / D


logml_truth = numeric(length(D_vec))
approx = matrix(0, N_REP, length(D_vec))

i = 1
for (D in D_vec) {
  
  J = D * KAPPA
  
  ## priors --------------------------------------------------------------------
  mu_0 = rep(0, D)      # prior mean for beta
  tau  = 1 / 4          # precision: inverse of variance
  sigmasq = 4           # true variance (1 x 1) 
  
  ## true beta -----------------------------------------------------------------
  beta = sample(-10:10, D, replace = T) 
  
  ## simulate regression data --------------------------------------------------
  
  X   = matrix(rnorm(N * D), N, D)                # (N x D) design matrix
  eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))    # (N x 1) errors vector
  y   = X %*% beta + eps                          # (N x 1) response vector
  
  data = list(X = X, y = y)
  
  ## compute posterior parameters ----------------------------------------------
  
  Q_beta = 1 / sigmasq * (t(X) %*% X + tau * diag(1, D))
  Q_beta_inv = solve(Q_beta)
  b = 1 / sigmasq * t(X) %*% y
  mu_beta = Q_beta_inv %*% b
  
  prior = list(y = y, X = X, sigmasq = sigmasq, tau = tau, N = N, D = D)
  post = list(Q_beta = Q_beta, Q_beta_inv = Q_beta_inv, mu_beta = mu_beta, b = b)
  
  logml_truth[i] = lil(prior, post)   # -272.1202
  
  for (p in 1:N_REP) {

      u_samps = rmvnorm(J, mean = c(mu_beta), sigma = Q_beta_inv) %>% data.frame
      u_df = preprocess(u_samps, D, prior) # J x (D + 1) -- stored row-wise

      approx[p, i] = logml(u_df, post)
  }

  cat("D = ", D, ", ", "J = ", J, ";", "\t",
      "LIL = ", round(logml_truth[i], 4),
      " (", round(mean(approx[,i]), 4), ")", '\n', sep = "")
  
  i = i + 1
}


# plot average error vs dimension 


delta = sweep(approx, 2, logml_truth)

x11()
df = data.frame(d = D_vec, 
                mae = colMeans(abs(delta)), 
                ae = colMeans(delta))

ggplot(df, aes(d, ae)) + geom_point(size = 2) + 
  geom_hline(yintercept = 0, col = 'red', size = 1) + 
  labs(y = "average error", title = "MVN (Kappa = 4)")






