



# lambda_k : (p x 1) grad_psi(m_k), m_k is representative point of partition A_k
# part_k   : (p x 2) intervals whose cartesian product defines the partition
Z_hat = function(lambda_k, part_k) {
    
    a_k = part_k[,1]  # lower bound of the p intervals
    b_k = part_k[,2]  # upper bound of the p intervals
    
    int_k = - 1 / lambda_k * exp(lambda_k * a_k) / exp(lambda_k * b_k) # p x 1
    
    return(int_k)
}


lambda_k = rep(1,2)
part_k = data.frame(matrix(c(-Inf, Inf, 7.1, Inf), 2, 2, byrow = T))


Z_hat(lambda_k, part_k)


