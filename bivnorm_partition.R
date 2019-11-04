
# bivnorm_partition.R ----------------------------------------------------------

# bivariate normal example 



library("MASS")
library(mixtools)  #for ellipse

N <- 200 # Number of random samples
set.seed(123)
# Target parameters for univariate normal distributions
rho <- -0.6
mu1 <- 1; s1 <- 2
mu2 <- 1; s2 <- 8

# Parameters for bivariate normal distribution
mu <- c(mu1,mu2) # Mean
sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2),
                2) # Covariance matrix

# Function to draw ellipse for bivariate normal data
ellipse_bvn <- function(bvn, alpha){
    Xbar <- apply(bvn,2,mean)
    S <- cov(bvn)
    ellipse(Xbar, S, alpha = alpha, col="red")
}

bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
colnames(bvn1) <- c("bvn1_X1","bvn1_X2")


# bivariate normal simulation --------------------------------------------------
# ------------------------------------------------------------------------------

library(ggplot2)

psi_biv = function(x_in, mu, Sigma) {
    # log-likelihood using the (posterior) samples
    return(dmvnorm(x_in, mu = mu, sigma = Sigma))
}

# covariance matrix: Sigma = (d / n) * Omega
# n >> d to make the distribution more concentrated

N = 500
d = 2
xi = d / 10
Omega = matrix(c(1, 0.7, 0.7, 1), nrow = 2)

# 
x_samp = mvrnorm(N, mu = c(1, 2), Sigma = xi * Omega) # stored row-wise




x_df = data.frame(x1 = x_samp[,1], x2 = x_samp[,2])

ggplot(x_df, aes(x1, x2)) + geom_point() + 
    geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0))


psi_u = psi_biv(x_in)
u_df = data.frame(x1 = x_in[,1], x2 = x_in[,2], psi_u = psi_u)


# (3) fit the decision tree

# (4) plot the partition over the parameter space: x1, x2

# (5) look at posterior distributions of mu, sigma_sq to see if there are more
#     partitions in areas of higher posterior probability

hist(x1)

hist(x2)
















