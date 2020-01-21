

source("C:/Users/ericc/mlike_approx/partition/partition.R")
source("C:/Users/ericc/mlike_approx/hybrid_approx.R")
setwd("C:/Users/ericc/mlike_approx/mvn")
source("mvn_helper.R")


D  = 5
N  = 5000
Sigma = D / N * diag(1, D)
Sigma_inv = solve(Sigma)

prior = list(Sigma = Sigma, Sigma_inv = Sigma_inv)


D / 2 * log(2 * pi) + 0.5 * log_det(Sigma)

set.seed(1)
J = 6000
N_approx = 10
u_samps = rmvnorm(J, mean = rep(0, D), sigma = Sigma) %>% data.frame 
u_df_full = preprocess(u_samps, D, prior)
approx_skew = approx_lil(N_approx, D, u_df_full, J / N_approx, prior)
mean(approx_skew)



N_vec_log = seq(6, 12, 0.05)        # sample size that is uniform over log scale
N_vec     = floor(exp(N_vec_log))   # sample size to use to generate data
logZ_0    = numeric(length(N_vec))
for (i in 1:length(N_vec)) {
    
    
    N = N_vec[i]
    Sigma = D / N * diag(1, D)
    Sigma_inv = solve(Sigma)
    
    logZ_0[i] = D / 2 * log(2 * pi) + 0.5 * log_det(Sigma)
    
    #print(paste("iter = ", i, 
    #            " -- Calculating LIL for D = ", D, ", N = ", N, sep = ''))
}


lil_df = data.frame(logZ_0 = logZ_0, logn = log(N_vec))

formula1 = y ~ x
ggplot(lil_df, aes(logn, logZ_0)) + geom_point() + 
    labs(title = "") + 
    geom_smooth(method = lm, se = T, formula = formula1) +
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ_0~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16)





