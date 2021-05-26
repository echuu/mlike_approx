

set.seed(1)
dim_vec = c(13, 15)
i = 1
while (i <= 2) {

  print(paste('starting sims for ', dim_vec[i], ' stacks', sep = ''))

  ###### initialize model
  n_G5 = dim_vec[i] # number of G_5 graphs we want to stack
  G = diag(1, n_G5) %x% G_5
  p = ncol(G)
  V = n * diag(1, p)

  ## try computing the normalizing constant of G_9 first as sanity check
  FREE_PARAMS_ALL = c(upper.tri(diag(1, p), diag = T) & G)

  edgeInd = G[upper.tri(G, diag = TRUE)] %>% as.logical

  ## construct A matrix so that we can compute k_i
  A = (upper.tri(diag(1, p), diag = F) & G) + 0

  k_i  = colSums(A) # see step 2, p. 329 of Atay
  nu_i = rowSums(A) # see step 2, p. 329 of Atay
  b_i = nu_i + k_i + 1

  set.seed(1)
  Omega_G = rgwish(1, G, b, V) # generate the true precision matrix
  P = chol(solve(V)) # upper cholesky factor; D^(-1) = TT'  in Atay paper

  # params = list(G = G, P = P, p = p, edgeInd = edgeInd,
  #               b = b, nu_i = nu_i, b_i = b_i)
  # N = 0
  S = matrix(0, p, p)
  D = sum(edgeInd) # number of free parameters / dimension of parameter space
  index_mat = matrix(0, p, p)
  index_mat[upper.tri(index_mat, diag = T)][edgeInd] = 1:D
  # index_mat[upper.tri(index_mat, diag = T)]
  t_ind = which(index_mat!=0,arr.ind = T)

  index_mat[lower.tri(index_mat)] = NA
  vbar = which(index_mat==0,arr.ind = T) # non-free elements
  n_nonfree = nrow(vbar)

  params = list(G = G, P = P, p = p, D = D, edgeInd = edgeInd,
                b = b, nu_i = nu_i, b_i = b_i,
                t_ind = t_ind, n_nonfree = n_nonfree, vbar = vbar)


  Z_5  = log(2^(0.5*p1*b + 7)) + I_G(0.5*(b-2)) + (-0.5 * p1 * b - 7) * log(n)
  Z = n_G5 * Z_5

  grad = function(u, params) { fast_grad(u, params) }
  hess = function(u, params) { fast_hess(u, params) }
  J = 200
  samps = samplegw(J, G, b, N, V, S, solve(P), FREE_PARAMS_ALL)
  u_samps = samps$Psi_free %>% data.frame
  # u_df = preprocess(u_samps, D, params)     # J x (D_u + 1)
  u_df = gwish_preprocess(u_samps, D, params_G5)     # J x (D_u + 1)
  u_star = gwish_globalMode(u_df, params, params_G5)


  ###### run simulations -------------------------------------------------------

  n_sims       = 100
  hyb          = numeric(n_sims)
  # hyb_old      = numeric(n_sims)
  gnorm_approx = numeric(n_sims)
  bridge       = numeric(n_sims)
  bridge_warp  = numeric(n_sims)

  j = 1
  J = 1000
  set.seed(1)
  while (j <= n_sims) {

    samps = samplegw(J, G, b, N, V, S, solve(P), FREE_PARAMS_ALL)
    u_samps = samps$Psi_free %>% data.frame

    ### hyb estimator ----------------------------------------------------------
    u_df = gwish_preprocess(u_samps, D, params_G5)     # J x (D_u + 1)
    logzhat = hybml_gwish(u_df, params_G5, psi = psi, grad = grad, hess = hess,
                          u_0 = u_star)$logz
    hyb[j] = logzhat # hybrid

    ### bridge estimator -------------------------------------------------------
    u_samp = as.matrix(u_samps)
    colnames(u_samp) = names(u_df)[1:D]
    lb = rep(-Inf, D)
    ub = rep(Inf, D)
    names(lb) <- names(ub) <- colnames(u_samp)
    bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
                                                   log_posterior = log_density,
                                                   data = params_G5,
                                                   lb = lb, ub = ub,
                                                   silent = TRUE)
    bridge[j] = bridge_result$logml

    bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
                                                   log_posterior = log_density,
                                                   data = params_G5,
                                                   lb = lb, ub = ub,
                                                   method = 'warp3',
                                                   silent = TRUE)
    bridge_warp[j] = bridge_result$logml

    ### gnorm estimator --------------------------------------------------------
    gnorm_approx[j] = gnorm(G, b, V, J)


    ### display some information about the running avg of the estimator + error
    print(paste('iter ', j, ': ',
                'hyb = ',     round(mean(hyb[hyb!=0]), 3),
                ' (error = ', round(mean((Z - hyb[hyb!=0])), 3), '), ',
                'bse = ',  round(mean(bridge[bridge!=0]), 3),
                ' (error = ', round(mean((Z - bridge[bridge!=0])), 3), '), ',
                'wbse = ',  round(mean((bridge_warp[bridge_warp!=0])), 3),
                ' (error = ', round(mean(Z - bridge_warp[bridge_warp!=0]), 3), '), ',
                sep = ''))

    j = j + 1
  }



  ###### compute results -------------------------------------------------------

  truth = Z
  approx = data.frame(truth, hyb = hyb, gnorm = gnorm_approx, bridge = bridge,
                      bridge_warp = bridge_warp)
  approx_long = reshape2::melt(approx, id.vars = 'truth')

  truth = Z
  res_tbl =
    data.frame(logz      = colMeans(approx),
               approx_sd = apply(approx, 2, sd),
               avg_error = colMeans(truth - approx),            # avg error
               mae       = colMeans(abs(truth - approx)),       # mean absolute error
               rmse      = sqrt(colMeans((truth - approx)^2)))  # root mean square error




  ###### save model results ----------------------------------------------------

  model_name = paste("gw_", 'p_', p, '_D_', D, '_stack_', n_G5, sep = '')
  results_obj = list(name = model_name,
                     p = p,
                     D = D,
                     truth = truth,
                     results = res_tbl,
                     nMCMC = J,
                     n = n,
                     delta = b,
                     all_approx = approx,
                     Omega_G = Omega_G)
  save_loc = 'C:/Users/ericc/Documents/sim_results/'
  saveRDS(results_obj, file = paste(save_loc, model_name, ".rds", sep = ''))
  print(paste('saving resuls... ',
              paste(save_loc, model_name, ".rds", sep = ''), sep = ''))
  i = i + 1
}


m15 = "gw_p_75_D_180_stack_15"
g15 = readRDS(paste(save_loc, m15, ".rds", sep = ''))

m17 = "gw_p_85_D_204_stack_17"
g17 = readRDS(paste(save_loc, m17, ".rds", sep = ''))
print(paste("p = ", g17$p, ', D = ', g17$D, sep = ''))
g17$results

print(paste("p = ", p, ', D = ', D, sep = ''))
print(res_tbl)



