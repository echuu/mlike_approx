


# install.packages("BayesTree")

library("BayesTree")
library('mvtnorm')   # multivariate normal density
library('MASS')      # mvnorm()
library('ggplot2')
library('tree')

f = function(x){
    10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma = 1.0  #y = f(x) + sigma*z , z~N(0,1)
n = 100      #number of observations
set.seed(99)
x=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
Ey = f(x)
y=Ey+sigma*rnorm(n)
lmFit = lm(y~.,data.frame(x,y)) #compare lm fit to BART later
##run BART
set.seed(99)
bartFit = bart(x,y,ndpost=200) # default is ndpost=1000, this is to run example fast.
plot(bartFit) # plot bart fit
##compare BART fit to linear matter and truth = Ey
fitmat = cbind(y,Ey,lmFit$fitted,bartFit$yhat.train.mean)
colnames(fitmat) = c('y','Ey','lm','bart')
print(cor(fitmat))


set.seed(123)

d = 2
N = 100
mu1 = c(1, 2)
mu2 = c(4, 6)
Sigma1 = d / N * diag(c(1, 1))
Sigma2 = d / N * diag(c(1, 1))

# mixture weights
pi1 = 0.2
pi2 = 0.8

J = 5000

n1 = sum(runif(J) <= pi1) # draws from 1st gaussian
n2 = J - n1               # draw sfrom 2nd gaussian


# define psi() function --------------------------------------------------------
psi = function(u, mu1, mu2, Sigma1, Sigma2, pi1, pi2) {
    
    # log-density evaluated over the rows of u -- (J x 1)
    loglik = log(pi1 * dmvnorm(u, mu1, Sigma1) + 
                     pi2 * dmvnorm(u, mu2, Sigma2))
    
    return(loglik)
    
} # end psi() function


u = rbind(mvrnorm(n1, mu1, Sigma1), mvrnorm(n2, mu2, Sigma2)) # (J x 2)

# check distribution
u_df = data.frame(u1 = u[,1], u2 = u[,2])
# ggplot(u_df, aes(u1, u2)) + geom_point() + theme_bw()

# evaluate psi(u)
psi_u = psi(u, mu1, mu2, Sigma1, Sigma2, pi1, pi2)           # (J x 1)
u_df = data.frame(u1 = u_df$u1, u2 = u_df$u2, psi_u = psi_u) # (J x 3)

bartFit = bart(u_df[,1:2], u_df[,3], ndpost = 200)
plot(bartFit) # plot bart fit







