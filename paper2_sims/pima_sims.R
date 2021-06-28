
## pima_sims.R
## Simulations to compute the hybrid-ep estimate in the context of Pima Indians.
## We rely on implementations from Friel/Wyse and Moulines for the competiting
## estimators, and we simply append our results to their results to compare
## our performance

## see RunPimaModified.R in:
## 'normalizingconstant/Code Jason Wyse/RunPimaModified.R'
## which has M1, M2 initialization and global variable declarations, both of
## which are required forr any of the code below to run




## **Note** everything below needs to be re-run when we switch from M1 to M2
## or vice versa, since the dimension of the parameter space will have changed,
## so the values of d, D both have changed, as well as u_star. we also do the
## following 'test run' once so we pre-determine the value of u_star

out = E2$metropolis.hastings(Its = 250000, BurnIn = 50000, fix = fix)
samps = data.frame(unname(out)[,1:d])
u_df = preprocess(samps, d, NULL)
u_star = globalMode(u_df, NULL)
hybridml::hybml(u_df, NULL, grad = grad, hess = hess, u_0 = u_star)$logz

n_sims       = 10
hyb_m1       = numeric(n_sims) ## store log evidence for model 1
hyb_m2       = numeric(n_sims) ## store log evidence for model 2
# bridge       = numeric(n_sims)= numeric(n_sims)

j = 1
set.seed(1)
while (j <= n_sims) {

  ## change E1 <-> E2 depending on the model
  out = E2$metropolis.hastings(Its = 250000, BurnIn = 50000, fix = fix)
  samps = data.frame(unname(out)[,1:d])
  u_df = preprocess(samps, d, NULL)

  ### hyb estimator ------------------------------------------------------------
  res = hybridml::hybml(u_df, NULL, grad = grad, hess = hess, u_0 = u_star)
  res$logz
  hyb_m2[j] = res$logz # hybrid


  ### display some information about the running avg of the estimator
  print(paste('iter ', j, ': ',
              'hyb = ',     round(mean(hyb_m2[hyb_m2 != 0]), 3),
              sep = ''))
  j = j + 1
}


## testing laplace
E2$log.laplace.evidence()

## generating laplace estimator using grad, hess
u_star
0.5*(d)*log2Pi - 0.5*log_det(hess(u_star, params)) - psi(u_star)



m1_logz = hyb # evidence of model 2
mean(m1_logz)

m2_logz = hyb_m2 # evidence of model 1
mean(m2_logz)

# compute Bayes factor
bf_21 = exp(m1_logz - m2_logz)

# mean Bayes factor
mean(bf_21)

colMeans()

logz_m1 = data.frame(cbind(tab.M1, HYB_EP = m1_logz))
row.names(logz_m1) = NULL
logz_m1
names(logz_m1)[7] = "HYB-EP"
names(logz_m1)[2] = "L-MAP"


logz_m2 = data.frame(cbind(tab.M2, HYB_EP = m2_logz))
row.names(logz_m2) = NULL
logz_m2
names(logz_m2)[7] = "HYB-EP"
names(logz_m2)[2] = "L-MAP"

# par(mfrow=c(2,1))
# boxplot(logz_m1)
# title("Model 1")
# boxplot(logz_m2)
# title("Model 2")

library(ggplot2)
m1_long = reshape2::melt(logz_m1)
m1_plot = ggplot(m1_long, aes(x = variable, y = value)) +
  geom_boxplot(lwd = 0.9) +
  labs(x = '', y = '', title = 'Model 1') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        axis.text=element_text(size=22))
m1_plot


m2_long = reshape2::melt(logz_m2)
m2_plot = ggplot(m2_long, aes(x = variable, y = value)) +
  geom_boxplot(lwd = 0.9) +
  labs(x = '', y = '', title = 'Model 2') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        axis.text=element_text(size=22))
m2_plot


multiplot(m1_plot, m2_plot, cols = 2)



## compute BF_12 of the table from moulines
BF_mat = exp(logz_m1 - logz_m2)
colMeans(BF_mat)






