// [[Rcpp::depends(RcppEigen)]]
#include "fast_psi.h"

using namespace Rcpp;

// [[Rcpp::export]]
float psi(NumericVector u, List prior) {
	NumericMatrix y = prior["y"];
	NumericMatrix X = prior["X"];

	 float d = u.size() - 1;
	

	NumericVector mu_beta = prior["mu_beta"];
	MAT_TYPE V_beta = prior["V_beta"];
	NumericMatrix V_beta_inv = prior["V_beta_inv"];
	float a = prior["a_0"];
	float b = prior["b_0"];
	float sigmasq = u[d];

	float loglik = 0;
	
	float sd = std::sqrt(sigmasq);
	for (unsigned int i	= 0; i < X.rows(); i++) {
		float prod = 0;
		for (unsigned int j = 0; j < d; j++) {
			prod += X(i, j) * u[j];
		}
		
		loglik += dnorm_log(y(i, 0), prod, sd);
	}

	NumericVector diff = u[Range(0, d-1)]  - mu_beta;

	float prod = 0;
	for (unsigned int i = 0; i < diff.length(); i++)
		prod += diff[i] * diff[i];

	float logprior = a * std::log(b) - (d / 2.0) * std::log(2.0*M_PI) -
					 0.5 * std::log(V_beta.determinant()) - std::lgamma(a) -
					 (a + d / 2.0 + 1.0) * std::log(sigmasq) -
					 1.0 / sigmasq * (b + 0.5 * prod);

	return (-loglik - logprior);
}

// [[Rcpp::export]]
float dnorm_log(float x, float mean, float sd) {
	return std::log((1.0/(sd * std::sqrt(2.0*M_PI))) * (std::exp(-std::pow(x-mean, 2.0)/ (2.0*std::pow(sd, 2)))));
}
