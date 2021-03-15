
setwd("C:/Users/ericc/normalizingconstant/Code Jason Wyse")
source("ModelEvidence.R")

RadiataPine = read.table("RadiataPine.txt",sep = " ",header = TRUE)
y = RadiataPine$y
n = length(y)
X1 = cbind(rep(1,n),RadiataPine$x-mean(RadiataPine$x))

X2 = cbind(rep(1,n),RadiataPine$z - mean(RadiataPine$z))



X = X1 # (42 x 2)
# X = X2

mu0 = c(3000,185)
Lambda0 = diag(1,2)
Lambda0[1,1] = 0.06
Lambda0[2,2] = 6.00

a_0 = 3
b_0 = 2*300^2

tX = t(X)
d = ncol(X)
tau0 = Lambda0
alpha = 2 * a_0
delta = 2 * b_0


# n = length(DataObs)
# y = DataObs
# X = as.matrix(DataCovariates)
# tX = t(X)
# d = ncol(X)
# mu0 = PriorMeanMu
# tau0 = PriorSignalMu
# alpha = 2*PriorShapePrecision
# delta = 2*PriorRatePrecision

# convenient constants
logPi = log(pi)
log2Pi = log(2*pi)
XTX = tX%*%X
XTy = tX%*%y
M = XTX + tau0
cholM = chol(M)
log.detM = 2*sum(log(diag(cholM)))
invM = chol2inv(cholM)
choltau0 = chol(tau0)
invtau0 = chol2inv(choltau0)
log.dettau0 = 2*(sum(log(diag(choltau0))))
P = diag(1,n) - X%*%invM%*%tX
beta0 = invM%*%(tX%*%y + tau0%*%mu0)
yTy = t(y)%*%y
c0 = yTy + t(mu0)%*%(tau0%*%mu0) - t(beta0)%*%M%*%beta0
c1 = t(mu0)%*%tau0%*%mu0

LIL = -0.5*n*logPi + 0.5*log.dettau0 - 0.5*log.detM + 0.5*alpha*log(delta) + 
    lgamma((n+alpha)/2) - lgamma(alpha/2) -0.5*(n+alpha)*log(c0+delta)


params = list(Q_0 = Lambda0, mu0 = mu0, alpha = alpha, beta = beta,
              d = d, n = n, 
              X = X, y = y)


# define psi, gradient, hessian
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions
source("C:/Users/ericc/mlike_approx/EP/hybml.R")
source("C:/Users/ericc/mlike_approx/radiata/radiata_helper.R")
sourceCpp("C:/Users/ericc/mlike_approx/radiata/radiata.cpp")

params = list(Q_0 = Lambda0, mu0 = mu0, alpha = alpha, delta = delta,
              d = d, n = n, M = M, 
              X = X, y = y, Xty = t(X) %*% y,
              tau0 = Lambda0, beta0 = beta0,
              ldtau0 = log.dettau0)

n.its   = 505000 # num iterations to run MCMC
burn.in = 101000 # num burn-in
fix = list(); fix$vars = rep(FALSE, d + 1); fix$values = numeric(d + 1);
u_samps = gibbs_radiata(n.its,burn.in,fix)
u_samps = u_samps %>% as.data.frame # 200000 x 3
u_samps %>% dim

D = d + 1
u_df_all = preprocess(u_samps, D, params)
row.names(u_df_all) = NULL
u_df_all %>% head
u_df = u_df_all
# u_df = u_df_all[sample(1:nrow(u_df_all), 5000),]
u_df %>% head

LIL = -0.5*n*logPi + 0.5*log.dettau0 - 0.5*log.detM + 0.5*alpha*log(delta) + 
    lgamma((n+alpha)/2) - lgamma(alpha/2) -0.5*(n+alpha)*log(c0+delta)

hml_approx = hml_const(1, D, u_df, 1000, params)
hml_approx$const_vec


u_star = theta_star()

# hyb_map
MAP_LOC = which(u_df$psi_u == min(u_df$psi_u))
u_0 = u_df[MAP_LOC,1:D] %>% unname() %>% unlist()
hybml(u_df, params, psi, grad, hess, u_0 = u_0)

# hyb
u_0 = theta_star()
hybml(u_df, params, psi, grad, hess, u_0 = u_0)








