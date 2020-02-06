

library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio

setwd("C:/Users/ericc/mlike_approx/singular/test")


## using the posterior parameters we calculated in the mvn-ig example,
## sample from the posterior w/ STAN to see if the parameters are
## what we expect


## run data generation from mvn_ig.R -------------------------------------------



# sample from posterior via stan -----------------------------------------------

J         = 500          # number of MC samples per approximation
N_approx  = 1           # number of approximations
burn_in   = 2000         # number of burn in draws
n_chains  = 4            # number of markov chains to run
stan_seed = 123          # seed

J_iter = 1 / n_chains * N_approx * J + burn_in 


post_dat = list(p = p,
                a_n = a_n, b_n = b_n, 
                mu_star = c(mu_star), V_star = V_star)

mvnig_fit = stan(file   = 'mvn_ig_sampler.stan', 
                 data   = post_dat,
                 iter   = J_iter,
                 warmup = burn_in,
                 chains = n_chains,
                 seed   = stan_seed,
                 refresh = 0) # should give us J * N_approx draws

# skip to bottom ---------------------------------------------------------------

u_samp = rstan::extract(mvnig_fit, pars = c("beta", "sigmasq"), permuted = TRUE)

u_beta = u_samp$beta %>% data.frame()
u_sigmasq = u_samp$sigmasq

dim(u_beta) # (J * N_approx) x p

u_post = u_samp$beta %>% data.frame() # (J x 2) : post sample stored row-wise

# ggplot(u_post, aes(X1, X2)) + geom_point()

# ------------------------------------------------------------------------------


# evalute psi(u), for each of the posterior samples u --------------------------

u_post = data.frame(beta_post = u_beta, sigmasq_post = u_sigmasq)

psi_u = apply(u_post, 1, psi_true_mvn, post = post) %>% unname() # (J x 1)


# prepare input for the approximation scheme -----------------------------------

# (1.2) construct u_df -- this will require some automation for colnames
u_df_names = character(D + 1)
for (d in 1:D) {
    u_df_names[d] = paste("u", d, sep = '')
}
u_df_names[D + 1] = "psi_u"

# populate u_df
u_df = cbind(u_post, psi_u) # (J * N_approx) x (D + 1)

# rename columns (needed since these are referenced explicitly in partition.R)
names(u_df) = u_df_names


# compute approximation --------------------------------------------------------

# pre-compute psi(u) so the algorithm can directly start with fitting the tree

## start back here after fitting the stan object

u_df = preprocess(mvnig_fit, D, post)


u_rpart = rpart(psi_u ~ ., u_df)



stan_approx = approx_lil_stan(N_approx, prior, post, D, u_df, J)

stan_approx

mean(stan_approx, na.rm = TRUE) # -106.1804
sd(stan_approx, na.rm = TRUE)  # 4.491305



