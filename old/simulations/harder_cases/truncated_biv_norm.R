
# truncated_biv_norm.R ---------------------------------------------------------

# uncomment lines below to install missing packages
# install.packages("tmvtnorm")
# install.packages("data.table")
# install.packages("dplyr")
# install.packages("ggplot2")


library("tmvtnorm")   # truncated normal
library("data.table") # data handling
library("dplyr")      # data handling
library("ggplot2")    # visualizations
library("tree")


set.seed(123)

mu = c(0.5, 0.5)
Sigma = matrix(c(1, 0.8, 0.8, 2), 2, 2)

# number of samples to draw
J = 10000

# sample from bivariate normal
X = rmvnorm(J, mu, Sigma)
plot(X)

# truncation points: x1 truncated to (-1, 0.5), x2 truncated to (-inf, 4)
a = c(-1, -Inf) # lower truncation
b = c(0.5, 4)   # upper truncation

# define psi() function --------------------------------------------------------
psi = function(u, mu, Sigma, a, b) {
    loglik = dtmvnorm(u, mu, Sigma, a, b, log = TRUE)
    return(loglik)
} # end psi() function


# sample J = 1000 using rejection sampling
X_r = rtmvnorm(J, mean = mu, sigma = Sigma, lower = a, upper = b, 
               algorithm = 'rejection') # 1000 x 2
plot(X_r)

# sample J = 1000 using gibbs sampling --> we use this to fit the tree
X_g = rtmvnorm(J, mean = mu, sigma = Sigma, lower = a, upper = b, 
               algorithm = 'gibbs')     # 1000 x 2
plot(X_g)


# evaluate psi(u)
psi_u = psi(X_g, mu, Sigma, a, b)                            # (J x 1)
u_df = data.frame(u1 = X_g[,1], u2 = X_g[,2], psi_u = psi_u) # (J x 3)

# fit decision tree
u_tree = tree(psi_u ~ u1 + u2, u_df)


# overlay partition on scatterplot of points drawn from true density
plot(u_df[,1], u_df[,2], pch = 20, cex = 0.8, col = "pink",
     xlab = 'u1', ylab = 'u2', main = 'J = 1e4')
partition.tree(u_tree, add = TRUE, cex = 0.0001, ordvars = c("u1", "u2"))


# overlay the partition on scatterplot of nontruncated density
plot(X[,1], X[,2], pch = 20, cex = 0.8, col = "pink",
     xlab = 'u1', ylab = 'u2', main = 'truncated biv. normal, J = 1e4')
points(u_df[,1], u_df[,2], pch = 20, cex = 0.8, col = "cyan")
partition.tree(u_tree, add = TRUE, cex = 0.0001, ordvars = c("u1", "u2"))


# nicer looking overlay of the truncated density on the non-truncated density
X_all = bind_rows(data.frame(X), data.frame(X_g))
X_all$samp = "trunc"
X_all$samp[1:J] = "biv-norm"

ggplot(X_all, aes(X1, X2)) + geom_point(aes(col = X_all$samp))
