

library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio

setwd("C:/Users/ericc/mlike_approx/singular/test")


## using the posterior parameters we calculated in the mvn-ig example,
## sample from the posterior w/ STAN to see if the parameters are
## what we expect

post_dat = list(p = p,
                a_n = a_n, b_n = b_n, 
                mu_star = c(mu_star), V_star = V_star)

mvnig_fit = stan(file   = 'mvn_ig_sampler.stan', 
                 data   = post_dat,
                 iter   = 100 * 1000) # should give us 200,000 draws

u_samp = rstan::extract(mvnig_fit, pars = c("beta", "sigmasq"), permuted = TRUE)

dim(u_samp) # (2 * iter) x 2

u_post = u_samp$beta %>% data.frame() # (J x 2) : post sample stored row-wise

names(u_post) = c("sigmasq")

ggplot(u_post, aes(sigmasq)) + geom_histogram()
ggplot(u_post, aes(X1, X2)) + geom_point()

u_beta = u_samp$beta %>% data.frame()
u_sigmasq = u_samp$sigmasq


N_approx = 100

stan_approx = approx_lil_stan(N_approx, prior, post, D, u_beta, u_sigmasq)

mean(stan_approx, na.rm = TRUE) # -106.1804
var(stan_approx, na.rm = TRUE)  # 4.491305



