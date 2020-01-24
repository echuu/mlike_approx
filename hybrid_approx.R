
## hybrid_approx.R

library(dplyr)
library(rpart)

## available functions in this file
## 
## (1) log_det()
## (2) preprocess()
## (3) approx_lil()
##
## -----------------------------------------------------------------------------



## log_det() -------------------------------------------------------------------
## input :
##          xmat   : matrix
## output : 
##          (1 x 1) log(det(xmat))
##
log_det = function(xmat) {
    return(c(determinant(xmat, logarithm = T)$modulus))
}



## preprocess() ----------------------------------------------------------------
## input :
##          post_samps : posterior samples from gamma(u), stored row-wise
##          D          : dimension of parameter
##          prior      : parameters to be passed into psi(), lambda()
## output : 
##          u_df       : dataframe w/ one more column than post_samps, contains
##                       psi(u) evalued for each posterior sample
##
preprocess = function(post_samps, D, prior) {
    
    psi_u = apply(post_samps, 1, psi, prior = prior) %>% unname() # (J x 1)
    
    # (1.2) name columns so that values can be extracted by partition.R
    u_df_names = character(D + 1)
    for (d in 1:D) {
        u_df_names[d] = paste("u", d, sep = '')
    }
    u_df_names[D + 1] = "psi_u"
    
    # populate u_df
    u_df = cbind(post_samps, psi_u) # J x (D + 1)
    names(u_df) = u_df_names
    
    
    return(u_df)
    
} # end of preprocess() function



## approx_lil() ----------------------------------------------------------------
## input :
##          N_approx   : # of approximations to form
##          D          : dimension of parameter
##          u_df_full  : (N_approx * J) x (D + 1) posterior samples, u,  and 
##                       psi(u) stored row-wise -- to be fed into rpart()
##          J          : # of MC samples to use PER APPROXIMATION
##          prior      : parameters to be passed into psi(), lambda()
## output : 
##          def_approx : (N_approx x 1) vector of approximations of LIL
##
approx_lil = function(N_approx, D, u_df_full, J, prior) {
    
    def_approx = numeric(N_approx)
    
    # compute approximation to LIL N_approx times
    for (t in 1:N_approx) {
        
        # if (t %% 10 == 0) {
        #     print(paste("iter", t))
        # }
        
        ## (1) subset out rows in u_df_full to be used in the t-th approximation
        row_id = J * (t - 1) + 1
        u_df = u_df_full[row_id:(row_id+J-1),]
        
        ## (2) fit the regression tree via rpart()
        u_rpart = rpart(psi_u ~ ., u_df)
        
        ## (3) process the fitted tree
        
        # (3.1) obtain the (data-defined) support for each of the parameters
        param_support = matrix(NA, D, 2) # store the parameter supports row-wise
        
        for (d in 1:D) {
            param_d_min = min(u_df[,d])
            param_d_max = max(u_df[,d])
            
            param_support[d,] = c(param_d_min, param_d_max)
        }
        
        # (3.2) obtain the partition
        u_partition = extractPartition(u_rpart, param_support)  # partition.R
        
        # organize all data into single data frame --> ready for approximation
        param_out = u_star(u_rpart, u_df, u_partition, D)
        
        n_partitions = nrow(u_partition)     # numebr of partitions 
        c_k          = numeric(n_partitions) # constant term for k-th partition
        zhat         = numeric(n_partitions) # integral over k-th partition
        
        # (4) compute closed form integral over each partition
        for (k in 1:n_partitions) {
            
            # extract "representative point" of the k-th partition
            star_ind = grep("_star", names(param_out))
            u = param_out[k, star_ind] %>% unlist %>% unname
            
            # compute lambda_k : gradient of psi, evaluated at u_star
            l_k = lambda(u, prior)       # (D x 1) 
            
            # evaluate e^c_k = e^{psi(u_star)}
            # c_k[k] = exp(-psi(u, prior)) # (1 x 1) - incorrect if NOT simple
            
            # compute e^{c_k}, the constant term in the integral for each k
            c_k[k] = exp(-psi(u, prior) + sum(l_k * u)) 

            # store each component of the D-dim integral 
            integral_d = numeric(D)      # (D x 1)
            
            for (d in 1:D) {
                
                # find column id of the first lower bound
                col_id_lb = grep("u1_lb", names(param_out)) + 2 * (d - 1)
                col_id_ub = col_id_lb + 1
                
                # d-th integral computed in closed form
                integral_d[d] = - 1 / l_k[d] * 
                    (exp(- l_k[d] * param_out[k, col_id_ub]) - 
                     exp(- l_k[d] * param_out[k, col_id_lb])) 
                
                # comment below to use the 2nd term in the taylor expansion
                # integral_d[d] = (param_out[k, col_id_ub] - 
                #                  param_out[k, col_id_lb])
                
            } # end of loop computing each of 1-dim integrals
            
            # compute the D-dim integral (product of D 1-dim integrals)
            zhat[k] = prod(c_k[k], integral_d)
            
        } # end of for loop over the K partitions
        
        # store the log integral \approx log marginal likelihood
        def_approx[t] = log(sum(zhat))
        
        # check for underflow, flip sign if underflow (this usually works)
        #if (is.nan(def_approx[t])) {
        #    def_approx[t] = log(-sum(zhat))
        #}
        
    } # end of N_approx outer loop
    
    return(def_approx)
    
} # end of approx_lil() function

