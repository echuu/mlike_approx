

setwd("C:/Users/ericc/normalizingconstant/Code Jason Wyse")
source("ModelEvidence.R")

RadiataPine = read.table("RadiataPine.txt",sep = " ",header = TRUE)
y = RadiataPine$y
n = length(y)
X1 = cbind(rep(1,n),RadiataPine$x-mean(RadiataPine$x))

X2 = cbind(rep(1,n),RadiataPine$z - mean(RadiataPine$z))



X = X1 # (42 x 2)
X = X2

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

hml_approx = hml_const(1, D, u_df, 1000, params)
hml_approx$const_vec


hybml(u_df, params, psi, grad, old_hess)

LIL

n_reps = 20
hyb_results = numeric(n_reps)
D = d + 1

for (i in 1:n_reps) {
    u_samps = gibbs_radiata(n.its,burn.in,fix) %>% as.data.frame
    u_df = preprocess(u_samps, D, params)
    row.names(u_df) = NULL
    hyb_results[i] = hybml(u_df, params, psi, grad, hess) %>% suppressWarnings()
    print(paste('iter: ', i, ' | hybrid = ', round(hyb_results[i], 3), 
                ' (error = ', round(hyb_results[i] - LIL, 3), ')',
                sep = ''))
}

x11()
results = data.frame(LIL = LIL, hyb_results = hyb_results, misc = hyb_results)
results_long = melt(results, id.vars = "LIL")
results_long %>% head
ggplot(results_long, aes(x = variable, y = value)) + geom_boxplot() + 
    geom_hline(yintercept = LIL, col = 'red', size = 1, linetype = 'dashed')
hyb_results


ggplot(delta_df, aes(x = variable, y = value)) + geom_boxplot() +
    geom_hline(yintercept = 0, col = 'red', size = 1, linetype = 'dashed') + 
    coord_flip() + 
    labs(y = expression(paste(Delta, ' ', ln, ' ', p(y))), x = '') +
    theme_bw() + 
    theme(axis.text  = element_text(size=25),
          axis.title = element_text(size=25,face="bold")) + 
    scale_y_continuous(breaks = seq(-60, 0, 10))

library(Rcpp)

params = list(Q_0 = Lambda0, mu0 = mu0, alpha = alpha, delta = delta,
              d = d, n = n, 
              X = X, y = y,
              tau0 = Lambda0, beta0 = beta0,
              ldtau0 = log.dettau0)

u1 = u_df[2,1:D] %>% unname %>% unlist

sourceCpp("C:/Users/ericc/mlike_approx/radiata/radiata.cpp")

grad(u1, params)
old_grad(u1, params)


microbenchmark(f1 = grad(u1, params),
               f2 = old_grad(u1, params))

out = grad(u1, params)
out %>% dim




