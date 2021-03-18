

library(TruncatedNormal)
library(tmg)
library(mvtnorm)

library(Rcpp)
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")     
source("C:/Users/ericc/mlike_approx/truncate/regTN_helper.R")

sourceCpp("C:/Users/ericc/mlike_approx/speedup/trunc_psi.cpp")



set.seed(123)
D = 100
N = 100
I_D = diag(1, D)

n_samps = 10
J       = 500
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

# write.csv(X,   "C:/Users/ericc/X.csv", row.names = FALSE)
# write.csv(y,   "C:/Users/ericc/y.csv", row.names = FALSE)
# write.csv(eps, "C:/Users/ericc/eps.csv", row.names = FALSE)
# write.csv(beta, "C:/Users/ericc/beta.csv", row.names = FALSE)


# beta_df = data.frame(beta, comp = 1:D)
# ggplot(beta_df, aes(comp, y = 1, fill = beta)) +
#     geom_tile(color = 'black', size = 1) +
#     #geom_text(aes(label = round(beta, 4)), size = 5) +
#     scale_fill_gradient(low = "white", high = "red",
#                         limits=c(0, 1), 
#                         breaks = c(0, 0.5, 1)) + coord_fixed() + 
#     theme_bw() + labs(y = '', x = expression(beta)) +
#     scale_x_continuous(breaks = 1:D, name = expression(beta)) + 
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           panel.background = element_blank(), axis.line = element_line(colour = "black"),
#           axis.text.y=element_blank(),
#           legend.title = element_blank(),
#           text = element_text(size = 15))


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

# write.csv(u_df, "C:/Users/ericc/u_df.csv", row.names = FALSE)

hybrid = hybrid_ml(D, u_df, J, prior)
hybrid$zhat

lb = rep(0, D)
ub = rep(Inf, D)
colnames(samples) = names(u_df)[1:D]
names(lb) <- names(ub) <- colnames(samples)

log_density = function(u, data) {
    -psi(u, data)
}

library(bridgesampling)
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


came_approx(u_df, prior, post, J, D)

B = 20
hyb = numeric(B)
bridge = numeric(B)
wbse = numeric(B)
hme = numeric(B)
came = numeric(B)
set.seed(1)
i = 1
J = 75

while (i <= B) {
    
    # sample from posterior
    samples = rtmvnorm(J, c(mu_beta), Q_beta_inv, rep(0, D), rep(Inf, D))
    u_df = preprocess(data.frame(samples), D, prior)
    
    #### (1) hybrid estimator
    hybrid = hybrid_ml(D, u_df, J, prior)
    if (any(is.na(hybrid))) {print(paste("error in iteration", i)); next;}
    hyb[i] = hybrid$zhat

    lb = rep(0, D)
    ub = rep(Inf, D)
    colnames(samples) = names(u_df)[1:D]
    names(lb) <- names(ub) <- colnames(samples)
    bridge_result <- bridgesampling::bridge_sampler(samples = samples,
                                                    log_posterior = log_density,
                                                    data = prior,
                                                    lb = lb, ub = ub,
                                                    silent = TRUE)
    bridge[i] = bridge_result$logml
    
    # bridge_result <- bridgesampling::bridge_sampler(samples = samples,
    #                                                 log_posterior = log_density,
    #                                                 data = prior,
    #                                                 lb = lb, ub = ub,
    #                                                 silent = TRUE,
    #                                                 method = 'warp3')
    # 
    # wbse[i] = bridge_result$logml
    # hme[i] = hme_approx(u_df, prior, J, D, N)
    came[i] = came_approx(u_df, prior, post, J, D)
    
    print(paste("iter ", i, ': ',
                "hybrid = ", round(mean(hyb[hyb!=0]), 3),
                '; ', "avg_err = ", round(mean((LIL - hyb[hyb!=0])), 4),
                '; ', "abs_err = ", round(mean(abs(LIL - hyb[hyb!=0])), 4),
                ' // ',
                "bridge = ", round(mean(bridge[bridge!=0]), 3), '; ',
                "avg_err = ", round(mean((LIL - bridge[bridge!=0])), 4), '; ',
                "abs_err = ", round(mean(abs(LIL - bridge[bridge!=0])), 4),
                # ' // ',
                # "wbse = ", round(mean(wbse[wbse!=0]), 3), '; ',
                # "mae = ", round(mean((LIL - wbse[wbse!=0])), 4),
                sep = '')) 
    
    i = i + 1
}

#### read in NSE stuff
# nse stuff
data_file_loc = 'C:/Users/ericc/Dropbox/eric chuu research/aistats/rdata_files/'

nse = read.csv(paste(data_file_loc, 'nse_trunc_d20.csv', sep = ''),
               header = F)[,1]
nse %>% head

approx = data.frame(LIL, hyb = hyb, bridge = bridge, wbse = wbse,
                    came = came, hme = hme, nse = nse)
error = data.frame(approx = colMeans(approx), 
                   approx_sd = apply(approx, 2, sd),
                   mae = colMeans(abs(LIL - approx)),
                   ae = colMeans((LIL - approx)),
                   rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)
error
approx %>% head
names(approx) = c("LIL", "HybE", "BSE", "WBSE", "CAME", "HME", "NSE")

# ------------------------------------------------------------------------------

fig_loc = 'C:/Users/ericc/mlike_approx/final_sims/'
saveRDS(list(J = J, D = D, N = N, approx = approx, error = error),
        file = paste(fig_loc, 'trunc_d20_j45.RData', sep = ''))


trunc = readRDS(paste(fig_loc, 'trunc_d20_j45.RData', sep = ''))



delta_df = (LIL - approx) %>% dplyr::select(-c('LIL')) %>% melt()

delta_df %>% head

### figure dimensions for saving: (954 x 488) 
ggplot(delta_df, aes(x = variable, y = value)) + geom_boxplot() +
    geom_hline(yintercept = 0, col = 'red', size = 1, linetype = 'dashed') +
    coord_flip() + 
    labs(y = expression(paste(Delta, ' ', ln, ' ', p(y))), x = '') +
    theme_bw() + 
    theme(axis.text  = element_text(size=25),
          axis.title = element_text(size=25,face="bold")) + 
    scale_y_continuous(breaks = seq(-30, 0, 10))




fig_loc = 'C:/Users/ericc/mlike_approx/final_sims/'
saveRDS(list(J = J, D = D, N = N, approx_df = approx, error = error),
        file = paste(fig_loc, 'trunc_d20_j1e5.RData', sep = ''))


approx = trunc$approx

LIL = approx$LIL[1]

approx %>% head
delta_df = (LIL - approx) %>% dplyr::select(-c('LIL', 'NSE')) %>% melt()

delta_df %>% head



ggplot(delta_df, aes(x = variable, y = value)) + 
    geom_boxplot(lwd = 1, fatten = 1.5) +
    geom_hline(yintercept = 0, col = 'red', size = 1.5, linetype = 'dashed') +
    coord_flip() + 
    labs(y = expression(paste(Delta, ' ', ln, ' ', p(y))), x = '') +
    theme_bw() + 
    theme(axis.text  = element_text(size=25),
          axis.title = element_text(size=25,face="bold")) + 
    scale_y_continuous(breaks = seq(-30, 0, 10))




### figure dimensions for saving: (954 x 488) 
ggplot(delta_df, aes(x = variable, y = value)) + geom_boxplot() +
    geom_hline(yintercept = 0, col = 'red', size = 1, linetype = 'dashed') +
    coord_flip() + 
    labs(y = expression(paste(Delta, ' ', ln, ' ', p(y))), x = '') +
    theme_bw() + 
    theme(axis.text  = element_text(size=25),
          axis.title = element_text(size=25,face="bold")) + 
    scale_y_continuous(breaks = seq(-30, 0, 10))

fig_loc = 'C:/Users/ericc/mlike_approx/final_sims/'
saveRDS(list(J = J, D = D, N = N, approx_df = approx, error = error),
        file = paste(fig_loc, 'trunc_d20_j1e5.RData', sep = ''))

test_read = readRDS(paste(fig_loc, 'trunc_d20_j1e5.RData', sep = ''))
test_read$J
test_read$error
















