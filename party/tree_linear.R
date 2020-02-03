


source("partition/partition.R")  # load partition extraction functions
source("extractPartition.R")
source("hybrid_approx.R")        # load main algorithm functions
source("mvn/mvn_helper.R")       # load psi(), lambda()


x11()

library(pracma)

fun = function(x, y) exp(-1/2 * (x^2 + y^2))

# DELL_PATH = "C:/Users/chuu/mlike_approx"
# setwd(DELL_PATH)

# model settings ---------------------------------------------------------------

D  = 2
Sigma = diag(1, D)
Sigma_inv = solve(Sigma)
prior = list(Sigma = Sigma, Sigma_inv = Sigma_inv)
D / 2 * log(2 * pi) + 0.5 * log_det(Sigma) # -3.683584, for D = 2, N = 500


set.seed(1)
J = 5000
N_approx = 1
u_samps = rmvnorm(J, mean = rep(0, D), sigma = Sigma) %>% data.frame 
u_df_full = preprocess(u_samps, D, prior)

# ------------------------------------------------------------------------------

x11()

# fit the tree using party library instead
# library(party)
install.packages("partykit", dependencies = TRUE)
library("partykit")

data("BostonHousing", package = "mlbench")
BostonHousing$lstat <- log(BostonHousing$lstat)
BostonHousing$rm <- BostonHousing$rm^2
BostonHousing$chas <- factor(BostonHousing$chas, levels = 0:1, 
                            labels = c("no", "yes"))
BostonHousing$rad <- factor(BostonHousing$rad, ordered = TRUE)

# these settings are all defaults exc for minsplit = 40 and verbose = TRUE
ctrl <- mob_control(alpha = 0.05, bonferroni = TRUE, minsplit = 40,
                    objfun = deviance, verbose = TRUE)

bhtree <- lmtree(medv ~ crim + zn | crim + zn, 
                 data = BostonHousing, minsize = 40)

print(bhtree)
plot(bhtree)

partykit:::.list.rules.party(bhtree)



u_tree = lmtree(psi_u ~ u1 | u2, data = u_df_full, minsize = 40)

partykit:::.list.rules.party(u_tree)

list.rules.party(u_tree)




















