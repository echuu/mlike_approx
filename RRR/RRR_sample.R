
##  RRR_sample.R 
##  read in the posterior samples written to .csv from MATLAB script

# path for lenovo
LEN_PATH  = "C:/Users/ericc/mlike_approx"
setwd(LEN_PATH)

# files must be loaded in this order, since *_helper.R file will sometimes
# overwrite a functino defined in hybrid_approx.R depending on the example

source("partition/partition.R")         # load partition extraction functions
source("hybrid_approx.R")               # load main algorithm functions
source("mvn/mvn_helper.R")              # load psi(), lambda() function
source('extractPartition.R')            # load extractPartition() function

# various dimension settings
p = 8            # number of columns in X
q = 6            # number of columns in Y
r = 2            # number of columns in B and A
n = 100          # number of rows in X and Y
sig2 = 10^(-2);  # fixed for now

# prior parameters
del = 10^(-2);


# dimension of the parameter, u (drawn from posterior)
d = r * p + r * q


u_df = read.csv("RRR/u_df_rrr.csv", header = FALSE)

dim(u_df)




# obtain A, B^T

# A_g  = reshape(u_df(g,1:p*r), p, r);    # recover the matrix A
# Bt_g = reshape(u_df(g,p*r+1:D), r, q);  # recover the matrix B^T =: B


# A is (p x r)
matrix(u_df[1000, 1:(p*r)], p, r) # fills the matrix COLUMN-wise (what we want)

# B^T is (r x q)
matrix(u_df[1000, (p*r+1):d], r, q)
