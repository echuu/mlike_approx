#include "fast_psi.h"

using namespace Rcpp;

// [[Rcpp::export]]
float psi(NumericVector u, List prior) {
    NumericMatrix y = prior["y"];
    NumericMatrix X = prior["X"];
    NumericVector mu_beta = prior["mu_beta"];
    
    float a = prior["a_0"];
    float b = prior["b_0"];
    float d = u.size() - 1;
    float sigmasq = u[d];
    float term1 = std::log((1.0/(std::sqrt(2.0 * M_PI * sigmasq))));
    float sigma2 = 2.0 * sigmasq;
    float loglik = 0;
    
    unsigned short x_lim = X.rows();	
    for (unsigned int i	= 0; i < x_lim; i++) {
        float prod = 0;
        for (unsigned int j = 0; j < d; j++) {
            prod += X(i, j) * u[j];
        }
        
        //loglik += dnorm_log(y(i, 0), prod, sigma2, term1);
        loglik += (term1 - std::pow(y(i, 0) - prod, 2.0) / sigma2);
    }
    
    NumericVector diff = u[Range(0, d-1)]  - mu_beta;
    
    float prod = 0;
    for (unsigned int i = 0; i < diff.length(); i++)
        prod += diff[i] * diff[i];
    
    float logprior = a * std::log(b) - (d / 2.0) * std::log(2.0 * M_PI) +
        std::lgamma(a) -
        (a + d / 2.0 + 1.0) * std::log(sigmasq) -
        1.0 / sigmasq * (b + 0.5 * prod);
    
    
    return (-loglik - logprior);
}

float dnorm_log(float x, float mean, float sigma2, float term1) {
    return term1 - std::pow(x-mean, 2.0) / sigma2;
}