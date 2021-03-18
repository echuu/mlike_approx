
library(bridgesampling)

sample_beta = function(s2, post) { 
    rmvnorm(1, mean = post$mu_star, sigma = s2 * post$V_star)
}

log_density = function(u, data) {
    -psi(u, data)
}


sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
beta_mat = t(sapply(sigmasq_post, sample_beta, post = post))
u_samp = data.frame(beta_mat, sigmasq_post)
u_df = preprocess(u_samp, D, prior)
hml_approx = hml_const(1, D, u_df, J, prior) 

lb <- c(rep(-Inf, p), 0)
ub <- c(rep(Inf, p), Inf)
colnames(u_samp) = names(u_df)[1:D]
names(lb) <- names(ub) <- colnames(u_samp)



B = 100
hme     = numeric(B) # store hybrid estimator
came    = numeric(B)
came0   = numeric(B)
bridge  = numeric(B) # bridge estimator (normal)

set.seed(1)
# i = 1
for (i in 1:B) {
    
    # sample from posterior
    sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
    beta_mat = t(sapply(sigmasq_post, sample_beta, post = post))
    u_samp = data.frame(beta_mat, sigmasq_post)
    u_df = preprocess(u_samp, D, prior)

    #### (1) harmonic mean estimator
    hme[i] = hme_approx(u_df, prior, J, D, N)
    
    ### (2) corrected arithmetic mean estimator (IS)
    hml_approx = hml_const(1, D, u_df, J, prior) 
    came_result = came_approx(u_df, hml_approx, prior, post, J, D)
    came[i] = came_result[1]
    came0[i] = came_result[2]

    #### (3) bridge sampling estimator
    u_samp = as.matrix(u_samp)
    colnames(u_samp) = names(u_df)[1:D]
    bridge_result = bridge_sampler(samples = u_samp,
                                   log_posterior = log_density,
                                   data = prior, lb = lb, ub = ub, 
                                   silent = TRUE)
    bridge[i] = bridge_result$logml
    
    print(paste("iter ", i, ': ',
                "bridge = ", round(mean(bridge[bridge!=0]), 3), '; ',
                "came = ",   round(mean(came0[came0!=0]), 3), '; ',
                "hme = ",    round(mean(hme[hme!=0]), 3), '; ',
                sep = ''))
}

LIL = lil(y, X, prior, post)
approx = data.frame(LIL, 
                    bridge = bridge[bridge!=0], 
                    came0  = came0[came0!=0],
                    came   = came[came!=0],
                    hme    = hme[hme!=0])
data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           ae = colMeans(LIL - approx),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)


#### merge with other results
setwd("C:/Users/ericc/Dropbox/eric chuu research/aistats/rdata_files")
mvnig = readRDS('mvnig_d20_n100.RData')
mvnig$approx_df

mvnig_all = mvnig$approx_df %>% 
    dplyr::mutate(bridge = bridge, came = came0, hme = hme)
mvnig_all %>% head

mvnig_df = mvnig_all %>% dplyr::mutate(iter = 1:nrow(mvnig_all)) %>% 
    dplyr::select(-c("LIL", "hme"))

mvnig_df_long = melt(mvnig_df, id.vars = 'iter')
mvnig_df_long %>% head # -303.8482

ggplot(mvnig_df, aes(x = iter, y = hyb)) + geom_point() +
    geom_hline(aes(yintercept = LIL), linetype = 'dashed', size = 0.9)

ggplot(mvnig_df_long, aes(x = iter, y = value, color = variable, shape = variable)) + 
    geom_point(size = 2) +
    geom_hline(aes(yintercept = LIL), linetype = 'dashed', size = 1.1) + 
    theme_bw() + 
    theme(legend.position = c(.9,.1)) + 
    theme(legend.title=element_blank()) + 
    scale_x_continuous(limits=c(0,100))
    
    
# image size: 780 x 512


save.image()






