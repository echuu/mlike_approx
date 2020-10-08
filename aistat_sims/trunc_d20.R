

library(TruncatedNormal)
library(tmg)
library(mvtnorm)

library(Rcpp)
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")     
source("C:/Users/ericc/mlike_approx/truncate/regTN_helper.R")

sourceCpp("C:/Users/ericc/mlike_approx/speedup/trunc_psi.cpp")



set.seed(123)
D = 20
N = 100
I_D = diag(1, D)

n_samps = 10
J       = 50
B       = 100 # number of replications

source("C:/Users/ericc/mlike_approx/paper_simulations/table2/mvn_estimators.R") 


#### generate data -------------------------------------------------------------

# prior mean
mu_0 = rep(0, D)
tau     = 1 / 4          # precision: inverse of variance
sigmasq = 4              # true variance (1 x 1) 

# true value of beta
set.seed(1)
beta = sample(0:10, D, replace = T)
beta = runif(D)
beta = c(runif(D-1, 0, 1), 0)

# generate the regression data -------------------------------------------------

X   = matrix(rnorm(N * D), N, D)                # (N x D) design matrix
eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))    # (N x 1) errors vector
y   = X %*% beta + eps                          # (N x 1) response vector

write.csv(X,   "C:/Users/ericc/X.csv", row.names = FALSE)
write.csv(y,   "C:/Users/ericc/y.csv", row.names = FALSE)
write.csv(eps, "C:/Users/ericc/eps.csv", row.names = FALSE)
write.csv(beta, "C:/Users/ericc/beta.csv", row.names = FALSE)



# compute posterior parameters -------------------------------------------------
Q_beta     =  1 / sigmasq * (t(X) %*% X + tau * diag(1, D))
Q_beta_inv =  solve(Q_beta)
b          =  1 / sigmasq * t(X) %*% y
mu_beta    =  Q_beta_inv %*% b


# create prior, post objects to be passed into the hml algorithm 
prior = list(y = y, X = X, sigmasq = sigmasq, tau = tau, N = N, D = D,
             Q_beta = Q_beta, b = b)
post = list(Q_beta = Q_beta, Q_beta_inv = Q_beta_inv, mu_beta = mu_beta, b = b)



#### compute true log marginal likelihood --------------------------------------

lil_0 = -0.5 * N * log(2 * pi) - 0.5 * (N + D) * log(sigmasq) + 
    0.5 * D * log(tau) - 0.5 * log_det(Q_beta) - 
    1 / (2 * sigmasq) * sum(y^2) + 0.5 * sum(b * mu_beta)

(LIL  = lil_0 + D * log(2) + 
        log(TruncatedNormal::pmvnorm(mu_beta, Q_beta_inv, 
                                     lb = rep(0, D), ub = rep(Inf, D))[1]))

samples = rtmvnorm(J, c(mu_beta), Q_beta_inv, rep(0, D), rep(Inf, D))
u_df = preprocess(data.frame(samples), D, prior)

write.csv(u_df, "C:/Users/ericc/u_df.csv", row.names = FALSE)

hybrid = hybrid_ml(D, u_df, J, prior)
hybrid$zhat

lb = rep(0, D)
ub = rep(Inf, D)
colnames(samples) = names(u_df)[1:D]
names(lb) <- names(ub) <- colnames(samples)

log_density = function(u, data) {
    -psi(u, data)
}

bridge_result <- bridge_sampler(samples = samples, log_posterior = log_density,
                                data = prior, lb = lb, ub = ub, silent = TRUE)
bridge_result$logml

hme_approx(u_df, prior, J, D, N)


# ------------------------------------------------------------------------------
# CAME function
# param_support = extractSupport(u_df, D) #
came_approx = function(u_df, prior, post, J, D) {
    
    A_beta = extractSupport(u_df, D) # data defined support (D x 2)
    post_mean = unlist(unname(colMeans(u_df[,1:D])))
    post_cov = cov(u_df[,1:D])
    # imp_samp = rmvnorm(J, mean = post_mean, sigma = post$Q_beta_inv) 
    
    imp_samp = rtmvnorm(J, c(post_mean), post_cov, rep(0, D), rep(Inf, D))
    
    u_df_imp = preprocess(data.frame(imp_samp), D, prior)
    
    log_s_theta = dtmvnorm(imp_samp, post_mean, post_cov, 
                           lb, ub, log = T, type = c("mc", "qmc"), B = 10000)
    
    include_d = rep(TRUE, J)
    for (d in 1:D) {
        include_d = include_d & 
            (imp_samp[,d] >= A_beta[d,1]) & (imp_samp[,d] <= A_beta[d,2])
    }
    
    -log(sum(include_d)) + log_sum_exp((-u_df_imp$psi_u - log_s_theta)[include_d])
}

# samples = rtmvnorm(J, c(mu_beta), Q_beta_inv, rep(0, D), rep(Inf, D))
# u_df = preprocess(data.frame(samples), D, prior)
# hml_approx = hml_const(1, D, u_df, J, prior)
# hml_approx$const_vec


B = 100
hyb = numeric(B)
bridge = numeric(B)
hme = numeric(B)
came = numeric(B)
came0 = numeric(B)
set.seed(1)

for (i in 1:B) {
    
    # sample from posterior
    samples = rtmvnorm(J, c(mu_beta), Q_beta_inv, rep(0, D), rep(Inf, D))
    u_df = preprocess(data.frame(samples), D, prior)
    
    # #### (1) hybrid estimator
    # hybrid = hybrid_ml(D, u_df, J, prior)
    # if (any(is.na(hybrid))) {print(paste("error in iteration", i)); next;}
    # hyb[i] = hybrid$zhat
    # 
    # lb = rep(0, D)
    # ub = rep(Inf, D)
    # colnames(samples) = names(u_df)[1:D]
    # names(lb) <- names(ub) <- colnames(samples)
    # 
    # # bridge_result <- bridge_approx(samples, log_density, prior, lb, ub)
    # # if (!is.na(bridge_result)) {
    # #     bridge[i] = bridge_result
    # # }
    # 
    # bridge_result <- bridgesampling::bridge_sampler(samples = samples, 
    #                                                 log_posterior = log_density,
    #                                                 data = prior, 
    #                                                 lb = lb, ub = ub, 
    #                                                 silent = TRUE)
    # bridge[i] = bridge_result$logml
    
    # hme[i] = hme_approx(u_df, prior, J, D, N)
    came[i] = came_approx(u_df, prior, post, J, D)
    
    print(paste("iter ", i, ': ',
                # "hybrid = ", round(mean(hyb[hyb!=0]), 3), 
                # '; ', "mae = ", round(mean(abs(LIL - hyb[hyb!=0])), 4),
                # ' // ',
                # "bridge = ", round(mean(bridge[bridge!=0]), 3), '; ', 
                # "mae = ", round(mean(abs(LIL - bridge[bridge!=0])), 4),
                "came = ", round(mean(came[came!=0]), 3), '; ', 
                "mae = ", round(mean(abs(LIL - came[came!=0])), 4),
                sep = '')) 
}

approx = data.frame(LIL, hyb = hyb[hyb!=0], bridge = bridge[bridge!=0])
approx = data.frame(LIL, came = came[came!=0])
data.frame(approx = colMeans(approx), 
           approx_sd = apply(approx, 2, sd),
           mae = colMeans(abs(LIL - approx)),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)


approx_long = (approx[,-1] - LIL) %>% melt()
ggplot(approx_long, aes(x = variable, y = value)) + geom_boxplot() +
    geom_hline(yintercept = 0) + 
    scale_y_continuous(breaks = seq(-4, 12, 0.5))


saveRDS(list(J = J, D = D, N = N, approx_df = approx), 
        file = 'trunc_d20.RData')
mvnig_d20 = readRDS('trunc_d20.RData')



approx = approx %>% mutate(iter = 1:B)
approx_long = approx %>% melt(id.vars = 'iter')
ggplot(approx_long, aes(x = iter, y = value, col = variable)) + geom_point() +
    geom_hline(aes(yintercept = LIL), linetype = 'dashed', size = 0.9) + 
    theme(legend.position = c(.2,.85))


# bridge stuff
approx = data.frame(LIL, bridge = bridge[bridge!=0])

trunc = readRDS('trunc_d20_n100.RData')
trunc$approx_df

trunc_df = trunc$approx_df %>% dplyr::mutate(iter = 1:nrow(mvnig$approx_df))
mvnig_df %>% head

ggplot(mvnig_df, aes(x = iter, y = hyb)) + geom_point() +
    geom_hline(aes(yintercept = LIL), linetype = 'dashed', size = 0.9) + 
    theme(legend.position = c(.2,.85))


# nse stuff
data_file_loc = 'C:/Users/ericc/Dropbox/eric chuu research/aistats/rdata_files/'

nse = read.csv(paste(data_file_loc, 'nse_trunc_d20.csv', sep = ''),
               header = F)[,1]
nse %>% head

approx = data.frame(LIL, nse = nse)
data.frame(approx = colMeans(approx), 
           approx_sd = apply(approx, 2, sd),
           mae = colMeans(abs(LIL - approx)),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)







