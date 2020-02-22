


library("mvtnorm")           # for draws from multivariate normal

# path for lenovo
LEN_PATH  = "C:/Users/ericc/mlike_approx"
setwd(LEN_PATH)

# path for dell
# DELL_PATH = "C:/Users/chuu/mlike_approx"
# setwd(DELL_PATH)

# files must be loaded in this order, since *_helper.R file will sometimes
# overwrite a functino defined in hybrid_approx.R depending on the example

source("partition/partition.R")         # load partition extraction functions
source("hybrid_approx.R")               # load main algorithm functions
source("mvn/mvn_helper.R")              # load psi(), lambda() function
source('extractPartition.R')            # load extractPartition() function

D  = 3
Omega = diag(1, D)
N  = 1000
Sigma = D / N * Omega 
Sigma_inv = solve(Sigma)
mu_0 = rep(0, D)
prior = list(Sigma = Sigma, Sigma_inv = Sigma_inv, mu_0 = mu_0)

# closed form of the normalizing constant
(logZ_mvn = D / 2 * log(2 * pi) + 0.5 * log_det(Sigma))


B        = 100          # number of batch estimates
N_approx = 1            # number of estimates to compute per iteration b


set.seed(1)

J       = 1000          # number of MC samples
u_samps = rmvnorm(J * B, mean = rep(0, D), sigma = Sigma) %>% data.frame 
u_df = preprocess(u_samps, D, prior)
u_df %>% head

hml_approx = hml(N_approx, D, u_df, J, prior)
hml_approx$hybrid_vec

hml_approx$n_const / nrow(hml_approx$verbose_partition)

# plotPartition(u_df, hml_approx$verbose_partition)

(logZ_mvn = D / 2 * log(2 * pi) + 0.5 * log_det(Sigma))



# ------------------------------------------------------------------------------


lil_hml  = numeric(B) # store the harmonic mean approximate for each batch
prop_const = numeric(B) # proportion of partitions that use constant approx 

for (b_i in 1:B) {
    
    if (b_i %% 20 == 0) print(paste("iter: ", b_i, "/", B, sep = ""))
    
    start = ((b_i - 1) * J + 1)
    end = start + J - 1
    
    u_df_b = u_df[start:end,]
    
    # (1) compute hybrid app
    hml_approx = hml(N_approx, D, u_df_b, J, prior) 
    lil_hml[b_i] = hml_approx$hybrid_vec
    
    prop_const[b_i] = hml_approx$n_const / nrow(hml_approx$verbose_partition)
    
} # end of simulation outer loop


## plot results
hme_df = data.frame(mcmc = 1:B, hml = lil_hml)

# mean average error (AE, true - estimated)
(MAE = round(mean(logZ_mvn - lil_hml), 3))

# root mean squared error (RMSE)
(RMSE = round(sqrt(mean((logZ_mvn - lil_hml)^2)), 3))

ggplot(hme_df, aes(x = mcmc, y = lil_hml)) + geom_point(col = 'blue') +
    geom_hline(aes(yintercept = logZ_mvn), linetype = 'dashed', size = 1.3) +
    labs(x = 'iter', y = 'logML',
         title = paste("logZ = ", round(logZ_mvn, 3), 
                       ", D = ", D, ", J = ", J, 
                       " (MAE = ", MAE, ", RMSE = ", RMSE, ")", sep = ''))




prop_const # prop of constant approximations used for each partition
mean(prop_const)




# ------------------------------------------------------------------------------

# perform analysis for a grid of J

J_vec   = c(100, 500, 1000, 2000)  # number of MC samples

lil_hml = data.frame(matrix(0, B, length(J_vec)))
names(lil_hml) = paste('j', J_vec, sep = '')

for (i in 1:length(J_vec)) {
    
    J = J_vec[i]
    
    u_samps = rmvnorm(J * B, mean = rep(0, D), sigma = Sigma) %>% data.frame 
    u_df = preprocess(u_samps, D, prior)
    
    for (b_i in 1:B) {
        
        if (b_i %% 20 == 0) print(paste("iter: ", b_i, "/", B, sep = ""))
        
        start = ((b_i - 1) * J + 1)
        end = start + J - 1
        
        u_df_b = u_df[start:end,]
        
        # (1) compute hybrid app
        hml_approx = hml(N_approx, D, u_df_b, J, prior) 
        lil_hml[b_i,i] = hml_approx$hybrid_vec
        
    } # end of simulation outer loop
    
}



## plot results
hme_df = cbind(mcmc = 1:B, lil_hml)

hme_df %>% head
# # mean average error (AE, true - estimated)
# (MAE = round(mean(logZ_mvn - lil_hml), 3))
# 
# # root mean squared error (RMSE)
# (RMSE = round(sqrt(mean((logZ_mvn - lil_hml)^2)), 3))

hme_df_long = melt(hme_df, id.vars = 'mcmc', value.name = 'logZ')


p1 = ggplot(hme_df, aes(x = mcmc, y = j100)) + geom_point(col = 'blue') +
    geom_hline(aes(yintercept = logZ_mvn), linetype = 'dashed', size = 1.3) +
    labs(x = 'iter', y = '',
         title = paste("logZ = ", round(logZ_mvn, 3), ", D = ", D, 
                       ", J = ", J_vec[1], 
                       ", approx = ", round(mean(hme_df$j100), 2), sep = '')) + 
    ylim(-6.5, -3.5) + 
    geom_hline(aes(yintercept = mean(hme_df$j100)), 
               col = 'red', linetype = 'dotdash', size = 1.3)
    

p2 = ggplot(hme_df, aes(x = mcmc, y = j500)) + geom_point(col = 'blue') +
    geom_hline(aes(yintercept = logZ_mvn), linetype = 'dashed', size = 1.3) +
    labs(x = 'iter', y = '',
         title = paste("logZ = ", round(logZ_mvn, 3), ", D = ", D, 
                       ", J = ", J_vec[2], 
                       ", approx = ", round(mean(hme_df$j500), 2), sep = '')) + 
    ylim(-6.5, -3.5) + 
    geom_hline(aes(yintercept = mean(hme_df$j500)), 
               col = 'red', linetype = 'dotdash', size = 1.3)


p3 = ggplot(hme_df, aes(x = mcmc, y = j1000)) + geom_point(col = 'blue') +
    geom_hline(aes(yintercept = logZ_mvn), linetype = 'dashed', size = 1.3) +
    labs(x = 'iter', y = '',
         title = paste("logZ = ", round(logZ_mvn, 3), ", D = ", D, 
                       ", J = ", J_vec[3], 
                       ", approx = ", round(mean(hme_df$j1000), 2), sep = '')) + 
    ylim(-6.5, -3.5) + 
    geom_hline(aes(yintercept = mean(hme_df$j1000)), 
               col = 'red', linetype = 'dotdash', size = 1.3)

p4 = ggplot(hme_df, aes(x = mcmc, y = j2000)) + geom_point(col = 'blue') +
    geom_hline(aes(yintercept = logZ_mvn), linetype = 'dashed', size = 1.3) +
    labs(x = 'iter', y = '',
         title = paste("logZ = ", round(logZ_mvn, 3), ", D = ", D, 
                       ", J = ", J_vec[4], 
                       ", approx = ", round(mean(hme_df$j2000), 2), sep = '')) + 
    ylim(-6.5, -3.5) + 
    geom_hline(aes(yintercept = mean(hme_df$j2000)), 
               col = 'red', linetype = 'dotdash', size = 1.3)

multiplot(p1, p2, p3, p4, cols = 4)






multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}





