
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")
source("C:/Users/ericc/mlike_approx/mvn_ig/mvn_ig_helper.R") 


library(Rcpp)
library(RcppEigen)
sourceCpp("C:/Users/ericc/mlike_approx/fast_psi.cpp")

# library(MLmetrics)

# load this LAST to overwrite def preprocess()


J = 100          # number of MC samples per approximation
D = 20
N = 100


set.seed(123)
p       = D - 1            # dimension of beta
mu_beta = rep(0, p)        # prior mean for beta
V_beta  = diag(1, p)       # scaled precision matrix for betaV_beta
a_0     = 2 / 2            # shape param for sigmasq
b_0     = 1 / 2            # scale param 
beta    = sample(-10:10, p, replace = T)
sigmasq = 4                # true variance (1 x 1) 

I_p = diag(1, p)           # (p x p) identity matrix
I_N = diag(1, N)           # (N x N) identity matrix

X = matrix(rnorm(N * p), N, p) # (N x p) design matrix

eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))

y = X %*% beta + eps # (N x 1) response vector
# ------------------------------------------------------------------

# write.csv(X, "C:/Users/ericc/X_mvnig.csv", row.names = F)
# write.csv(y, "C:/Users/ericc/y_mvnig.csv", row.names = F)
# 
# 
# 
# ## beta heat map
# beta_df = data.frame(beta, comp = 1:p)
# names(beta_df) = c(expression(beta), 'comp')
# ggplot(beta_df, aes(comp, y = 1, fill = beta)) +
#     geom_tile(color = 'black', size = 1) +
#     scale_fill_gradient(low = "white", high = "red",
#                         limits=c(-10, 10), 
#                         breaks = c(-10, -5, 0, 5, 10)) + coord_fixed() + 
#     theme_bw() + labs(y = '', x = expression(beta)) +
#     scale_x_continuous(breaks = 1:p, name = expression(beta)) + 
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           panel.background = element_blank(), axis.line = element_line(colour = "black"),
#           axis.text.y=element_blank(),
#           legend.title = element_blank(),
#           text = element_text(size = 15))

## compute posterior parameters ------------------------------------
V_beta_inv = solve(V_beta)
V_star_inv = t(X) %*% X + V_beta_inv

V_star  = solve(V_star_inv)                                # (p x p)
mu_star = V_star %*% (t(X) %*% y + V_beta_inv %*% mu_beta) # (p x 1)
a_n =  a_0 + N / 2 
b_n =  c(b_0 + 0.5 * (t(y) %*% y + 
                          t(mu_beta) %*% V_beta_inv %*% mu_beta - 
                          t(mu_star) %*% V_star_inv %*% mu_star))

# compute MLE estimates for mean, variance of the regression model
ybar = X %*% mu_star
sigmasq_mle = 1 / N * sum((y - ybar)^2)

# create prior, posterior objects
prior = list(V_beta = V_beta, mu_beta = mu_beta,  a_0 = a_0, b_0 = b_0,
             y = y, X = X, V_beta_inv = V_beta_inv)
# store posterior parameters
post  = list(V_star =  V_star, mu_star = mu_star, a_n = a_n, b_n = b_n, 
             V_star_inv = V_star_inv)

sample_beta = function(s2, post) { 
    rmvnorm(1, mean = post$mu_star, sigma = s2 * post$V_star)
}

(LIL = lil(y, X, prior, post))    # -256.7659


#### BRIDGE --------------------------------------------------------------------
library(bridgesampling)
log_density = function(u, data) {
    -psi(u, data)
}





#### CORRECTED ARITHMETIC MEAN ESTIMATOR ---------------------------------------
came_approx = function(u_df, prior, post, J, D) {
    
    A_supp = extractSupport(u_df, D) # data defined support (D x 2)
    
    ## draw from importance distribution
    post_mean = unlist(unname(colMeans(u_df[,1:D])))
    
    beta_mean    = post_mean[1:p]
    sigmasq_mean = post_mean[D]
    sigmasq_var  = var(u_df[,D])
    r_imp = sigmasq_mean^2/sigmasq_var + 2
    s_imp = sigmasq_mean * (1 + sigmasq_mean^2 / sigmasq_var)
    
    
    beta_imp = rmvnorm(J, mean = beta_mean, 
                       sigma = sigmasq_mean * post$V_star) %>% data.frame 
    sigmasq_imp = MCMCpack::rinvgamma(J, shape = r_imp, scale = s_imp)
    
    imp_samp = data.frame(beta_imp, sigmasq_imp)
    u_df_imp = preprocess(imp_samp, D, prior)
    
    # to be updated
    log_s_theta = unname(dmvnorm(beta_imp, mean = beta_mean, 
                                 sigma = sigmasq_mean * post$V_star, 
                                 log = TRUE)) + 
        log(MCMCpack::dinvgamma(sigmasq_imp, shape = r_imp, scale = s_imp))
    
    include_d = rep(TRUE, J)
    for (d in 1:D) {
        include_d = include_d & 
            (imp_samp[,d] >= A_supp[d,1]) & (imp_samp[,d] <= A_supp[d,2])
    }
    
    Jp = sum(include_d)
    
    # scame_j = -log(J) + log_sum_exp((-u_df_imp$psi_u - log_s_theta)[include_d])
    out = -log(Jp) + log_sum_exp((-u_df_imp$psi_u - log_s_theta)[include_d])
    
    return(out)
}





#### HARMONIC MEAN ESTIMATOR ---------------------------------------------------

# reg_lik() function -----------------------------------------------------------
# reformulate the likelihood so that it is of the form exp(x) so that 
# we can take advantage of the log-sum-exp trick (below); this function
# returns the part of the reformulated likelihood that is in the exponential
reg_lik = function(u_df, prior, J, D, N) {
    
    # lik_inv = numeric(J) # store likelihood computed for MCMC samples
    tmp_j   = numeric(J) # store quantity that is passed into log(sum(exp(x)))
    p = D - 1
    
    # sigmasq_samp = sigmasq
    
    for (j in 1:J) {
        
        beta_samp    = unname(unlist(u_df[j, 1:p]))
        # uncomment line below for NIG case
        sigmasq_samp = unname(unlist(u_df[j, D]))
        
        # lik_inv[j] = 1 / prod(dnorm(data$y, data$X %*% beta_samp, 
        #                         sqrt(sigmasq_samp)))
        
        tmp_j[j] = 1 / (2 * sigmasq_samp) * 
            sum((prior$y - prior$X %*% beta_samp)^2) + N / 2 * log(sigmasq_samp)
        
    } # end of loop iterating over MCMC samples
    
    
    return(tmp_j)
    
} # end of reg_lik() function --------------------------------------------------



# hme() function ---------------------------------------------------------------
# harmonic mean estimator -- this is written specifically for the MVN-IG
# example, since the likelihood is of a form such that we can take advantage
# of the log-sum-exp trick to stabilize the calculation of estimator
hme_approx = function(u_df, prior, J, D, N) {
    
    # harmonic mean estimator requires calculating the likelihood given
    # each of the J parameters (that are sampled via MCMC)
    
    # in order to generalize this function, each model that wants to take
    # advantage of the hme estimator should provide
    # (1) likelihood function
    # (2) parameter extraction that only requires u_df input
    
    # lik_inv = reg_lik(u_df, data, J, D)
    # hme_estimate = log(J) - log(sum(lik_inv))
    
    # log_sum_exp version of 
    tmp_j = reg_lik(u_df, prior, J, D, N)
    hme_estimate = log(J) - N / 2 * log(2 * pi) - log_sum_exp(tmp_j)
    
    return(hme_estimate)
    
} # end of hme() function ------------------------------------------------------



#### TEST ALL ESTIMATORS

# sample from posterior
sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
beta_mat = t(sapply(sigmasq_post, sample_beta, post = post))
u_samp = data.frame(beta_mat, sigmasq_post)
u_df = preprocess(u_samp, D, prior)

lb <- c(rep(-Inf, p), 0)
ub <- c(rep(Inf, p), Inf)
u_samp = as.matrix(u_samp)
colnames(u_samp) = names(u_df)[1:D]
names(lb) <- names(ub) <- colnames(u_samp)
bridge_result = bridge_sampler(samples = u_samp,
                               log_posterior = log_density,
                               data = prior, lb = lb, ub = ub, silent = TRUE)

(LIL = lil(y, X, prior, post))
bridge_result$logml
hybrid_ml(D, u_df, J, prior)$zhat
hme_approx(u_df, prior, J, D, N)
came_approx(u_df, prior, post, J, D)


B = 100
hyb    = numeric(B)
bridge = numeric(B)
hme    = numeric(B)
came   = numeric(B)
wbse   = numeric(B)
set.seed(1)
J = 100
i = 1

#### START SIMULATION LOOP -----------------------------------------------------
while (i <= B) {
    
    # sample from posterior
    sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
    beta_mat = t(sapply(sigmasq_post, sample_beta, post = post))
    u_samp = data.frame(beta_mat, sigmasq_post)
    u_df = preprocess(u_samp, D, prior)
    
    lb <- c(rep(-Inf, p), 0)
    ub <- c(rep(Inf, p), Inf)
    u_samp = as.matrix(u_samp)
    colnames(u_samp) = names(u_df)[1:D]
    names(lb) <- names(ub) <- colnames(u_samp)
    bridge_result = bridge_sampler(samples = u_samp,
                                   log_posterior = log_density,
                                   data = prior, lb = lb, ub = ub, silent = TRUE,
                                   method = 'warp3')
    # bridge[i] = bridge_result$logml
    wbse[i] = bridge_result$logml
    
    #### (1) hybrid estimator
    # hybrid = hybrid_ml(D, u_df, J, prior)
    # 
    # if (any(is.na(hybrid))) {print(paste("error in iteration", i)); next;}
    # 
    # hyb[i] = hybrid$zhat
    # 
    # came[i] = came_approx(u_df, prior, post, J, D)
    # hme[i] = hme_approx(u_df, prior, J, D, N)
    
    print(paste("iter ", i, ': ',
                "hybrid = ", round(mean(hyb[hyb!=0]), 3),
                '; ', "ae = ", round(mean((LIL - hyb[hyb!=0])), 4),
                ' // ',
                "bridge = ", round(mean(wbse[wbse!=0]), 3), '; ',
                "ae = ", round(mean((LIL - wbse[wbse!=0])), 4), '; ',
                "abs_error = ", round(mean(abs(LIL - wbse[wbse!=0])), 4),
                # "came = ", round(mean(came[came!=0]), 3), '; ', 
                # "ae = ", round(mean((LIL - came[came!=0])), 4),
                sep = '')) 
    i = i + 1
}

approx = data.frame(LIL, hyb = hyb[hyb!=0], bridge = bridge[bridge!=0], 
                    wbse = wbse,
                    came = came[came!=0], hme = hme[hme!=0])
(error = data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
                    ae = colMeans((LIL - approx)),
                    rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3))

fig_loc = 'C:/Users/ericc/mlike_approx/final_sims/'
saveRDS(list(J = J, D = D, N = N, approx = approx, error = error),
        file = paste(fig_loc, 'mvnig_d20_j45.RData', sep = ''))


colMeans((mvnig$approx_df - LIL))

LIL = mvnig$approx_df$LIL[1]

approx_df = mvnig$approx_df

names(approx) = c("LIL", "HybE", "BSE", "CAME", "HME")
delta_df = (LIL - approx) %>% dplyr::select(-c('LIL')) %>% melt()



delta_df %>% head

ggplot(delta_df, aes(x = variable, y = value)) + geom_boxplot(lwd = 0.8) +
    geom_hline(yintercept = 0, col = 'red', size = 1.5, linetype = 'dashed') +
    coord_flip() + 
    labs(y = expression(paste(Delta, ' ', ln, ' ', p(y))), x = '') +
    theme_bw() + 
    theme(axis.text  = element_text(size=25),
          axis.title = element_text(size=25,face="bold")) + 
    scale_y_continuous(breaks = seq(-30, 0, 10))






### figure dimensions for saving: (954 x 488) 
ggplot(delta_df, aes(x = variable, y = value)) + geom_boxplot() +
    coord_flip() + 
    labs(y = expression(paste(Delta, ' ', ln, ' ', p(y))), x = '') +
    theme_bw() + 
    theme(axis.text  = element_text(size=25),
          axis.title = element_text(size=25,face="bold")) + 
    scale_y_continuous(breaks = seq(-50, 0, 10))

approx = approx_df
(error = data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
                    ae = colMeans((LIL - approx)),
                    rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3))









