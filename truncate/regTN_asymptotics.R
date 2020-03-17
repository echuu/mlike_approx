







# ------------------------------------------------------------------------------

set.seed(1)
K_sims    = 10                               # num of sims to run for each N
N_vec_log = seq(5, 10, by = 0.25)            # sample size grid unif in log
N_vec     = floor(exp(N_vec_log)) %>% unique # sample size to generate data


for (i in 1:length(N_vec)) {
    
    N = N_vec[i]
    
    for (k in 1:K_sims) {
        
        # generate data --------------------------------------------------------
        
        
        # compute true LIL first
        # TODO: without forming approximation, plot the true LIL to verify
        # -D/2 slope
        
        
        # sample from posterior, u
        
        
        # evaluate psi(u)
        
        
        # compute approximation
        
        
        # subtract the loglikelihood
        
    }
    
    
    print(paste("approximating LIL for N = ", N, " -- LIL = ", 
                round(mean(LIL_N_k[i, k]), 2), sep = ''))
    
    
}

















library(reshape2)
library(ggpmisc)
formula1 = y ~ x

ggplot(LIL_df, aes(x = log_N, y = LIL_N)) + geom_point(size = 1.5) + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    labs(x = "log(n)", y = "log(ML)", 
         title = "Approximate ML for Logistic Regression, D = 3") + 
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16) + 
    theme(legend.position = "none")
