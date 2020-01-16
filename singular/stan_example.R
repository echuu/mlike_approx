

library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio

# visit url below for guide
# https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

# stan configuration 

dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, ifelse(.Platform$OS.type == "windows", "Makevars.win", "Makevars"))
if (!file.exists(M)) file.create(M)
cat("\nCXX14FLAGS=-O3 -march=native -mtune=native",
    if( grepl("^darwin", R.version$os)) "CXX14FLAGS += -arch x86_64 -ftemplate-depth-256" else 
        if (.Platform$OS.type == "windows") "CXX11FLAGS=-O3 -march=corei7 -mtune=corei7" else
            "CXX14FLAGS += -fPIC",
    file = M, sep = "\n", append = TRUE)

# enable option to run models in parallel
options(mc.cores = parallel::detectCores())

# save compiled stan program to hard disk so that don't need to recompile
rstan_options(auto_write = TRUE)

Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7')

# ------------------------------------------------------------------------------

setwd("C:/Users/chuu/mlike_approx/singular")

# school example

# the items of this list match the 'constructor' in the input data definition
# in 8schools.stan
school_dat = list(J = 8, 
                y = c(28,  8, -3,  7, -1,  1, 18, 12),
                sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

# compile 
fit = stan(file = '8schools.stan', data = school_dat)

# summary for the parameters of the model
# log posterior
print(fit)

# extract() function allows us to obtain samples 

# extract posterior samples for mu, tau
x = extract(fit, pars = c("mu", "tau"), permuted = TRUE)

head(x$mu)
head(x$tau)

# ------------------------------------------------------------------------------

N = 100
y = rnorm(N, mean = 13, sd = 2)

gamma_dat = list(N = N,
                 y = y)

gamma_fit = stan(file = 'gamma_sample.stan', data = gamma_dat)

u_samp = extract(gamma_fit, pars = c('mu', 'sigma'), permuted = TRUE)

# posterior samples of mu
head(u_samp$mu)

# posterior samples of sigma
head(u_samp$sigma)


## question: are mu, sigma above posterior samples? what are the priors?

# ------------------------------------------------------------------------------

# https://mc-stan.org/docs/2_21/reference-manual/increment-log-prob-section.html


N = 50
J = 2000


gamma_dat = list(N = N)

gamma_fit = stan(file = 'gamma_sample.stan', data = gamma_dat)

u_samp = extract(gamma_fit, pars = c("u"), permuted = TRUE)

head(u_samp$u %>% matrix())

u_test = u_samp$u[1:J,] %>% unname()

u1 = u_samp$u[1,]

