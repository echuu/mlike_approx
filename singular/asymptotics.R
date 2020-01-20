

#### does asymptotic analysis for a SINGLE value of D - see 
#### mvn_ig_asymptotics.R for analysis over a grid of D



library(mvtnorm)           # for draws from multivariate normal
library("numDeriv")        # for grad() function - numerical differentiation
library('MCMCpack')        # for rinvgamma() function
library('microbenchmark')

# path for lenovo
setwd("C:/Users/ericc/mlike_approx")

# path for dell
# setwd("C:/Users/chuu/mlike_approx")
source("partition/partition.R")      # load partition extraction functions
source("mvn_ig_helper.R")  # load functions specific to this model
setwd("C:/Users/ericc/mlike_approx/singular/test")


# STAN SETTINGS ----------------------------------------------------------------
J         = 100          # number of MC samples per approximation
N_approx  = 20           # number of approximations
burn_in   = 2000         # number of burn in draws, must be > (N_approx * J)
n_chains  = 4            # number of markov chains to run
stan_seed = 123          # seed

J_iter = 1 / n_chains * N_approx * J + burn_in 
# ------------------------------------------------------------------------------










