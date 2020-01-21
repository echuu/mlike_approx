

# psi() : negative log posterior
psi = function(u, prior) {
    
    N = prior$N
    
    return(N * u[1]^2 * u[2]^4)
}


# lambda() : gradient of psi
lambda = function(u, prior) {
    
    N = prior$N
    
    l_1 = 2 * N * u[1] * u[2]^4
    l_2 = 4 * N * u[1]^2 * u[2]^3
    
    return(c(l_1, l_2))
}



# 
# preprocess = function(post_samps, D, prior) {
#     
# 
#     psi_u = apply(post_samps, 1, psi, prior = prior) %>% unname() # (J x 1)
#     
#     # (1.2) construct u_df -- this will require some automation for colnames
#     u_df_names = character(D + 1)
#     for (d in 1:D) {
#         u_df_names[d] = paste("u", d, sep = '')
#     }
#     u_df_names[D + 1] = "psi_u"
#     
#     # populate u_df
#     u_df = cbind(u_post, psi_u) # (J * N_approx) x (D + 1)
#     
#     # rename columns (needed since these are referenced explicitly in partition.R)
#     names(u_df) = u_df_names
#     
#     
#     return(u_df)
#     
# }
# 
# 
# approx_lil_stan = function(N_approx, D, N, u_df_full, J, prior) {
#     
#     #### algorithm: main loop
#     # N_iters = N_approx
#     
#     # test_out = numeric()
#     def_approx = numeric(N_approx)  
#     
#     for (t in 1:N_approx) {
#         
#         
#         row_id = J * (t - 1) + 1
#         
#         u_df = u_df_full[row_id:(row_id+J-1),]
#         
#         ## (2) fit the regression tree via rpart()
#         u_rpart = rpart(psi_u ~ ., u_df)
#         # plot(u_rpart)
#         
#         ## (3) process the fitted tree
#         
#         # (3.1) obtain the (data-defined) support for each of the parameters
#         param_support = matrix(NA, D, 2) # store the parameter supports row-wise
#         
#         for (d in 1:D) {
#             # param_d_min = min(u_df[,d])
#             # param_d_max = max(u_df[,d])
#             param_d_min = 0
#             param_d_max = 1
#             
#             param_support[d,] = c(param_d_min, param_d_max)
#         }
#         
#         # paste code back here
#         # (3.2) obtain the partition --- moment of truth!!
#         u_partition = paramPartition(u_rpart, param_support)  # partition.R
#         
#         # organize all data into single data frame --> ready for approximation
#         param_out = u_star(u_rpart, u_df, u_partition, D)
#         
#         n_partitions = nrow(u_partition)
#         c_k = numeric(n_partitions)
#         zhat = numeric(n_partitions)
#         
#         for (k in 1:n_partitions) {
#             
#             star_ind = grep("_star", names(param_out))
#             u = param_out[k, star_ind] %>% unlist %>% unname
#             
#             c_k[k] = exp(-psi(u, prior)) # (1 x 1)
#             
#             l_k = lambda(u, prior)
#             
#             integral_d = numeric(D) # store each component of the D-dim integral 
# 
#             for (d in 1:D) {
# 
#                 # updated 1/14: find column id of the first lower bound
#                 col_id_lb = grep("u1_lb", names(param_out)) + 2 * (d - 1)
#                 col_id_ub = col_id_lb + 1
#                 
#                 # d-th integral computed in closed form
#                 integral_d[d] = - 1 / l_k[d] * 
#                     exp(- l_k[d] * (param_out[k, col_id_ub] - 
#                                         param_out[k, col_id_lb]))        
#                 
#             }
#             
#             zhat[k] = prod(c_k[k], integral_d)
#             
#         }
#         
#         def_approx[t] = log(sum(zhat))
#         
#         # if (is.nan(def_approx[t])) {
#         #     def_approx[t] = log(-sum(zhat))
#         # }
#         
#         
#         
#     }
#     
#     return(def_approx)
#     
#     # return(0)
#     
# } # end of approx_lil_stan()







