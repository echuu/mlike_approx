



logistic_lik = function(u_df, data, J, D, N) {
    
    lik = numeric(J)
    
    X = data$X
    y = data$y
    
    for (j in 1:J) {
    
        beta_samp = unname(unlist(u_df[j,1:D]))
        
        # lik[j] = prod(inv.logit(X %*% beta_samp)^y *
        #     (1 + exp(X %*% beta_samp))^(y - 1))
        
        lik[j] = -(sum(y * (X %*% beta_samp) - log(1 + exp(X %*% beta_samp))))
        
    }
    return(lik)
}


logistic_hme = function(u_df, data, J, D, N) {
    
    tmp_j = logistic_lik(u_df, data, J, D, N)
    
    hme_estimate = log(J) - log_sum_exp(tmp_j)
    
    return(hme_estimate)   
    
}



data = list(X = X, y = y)

lik = logistic_lik(u_df, data, J, D, N)

log(J) - log_sum_exp(lik)

logistic_hme(u_df, data, J, D, N)


