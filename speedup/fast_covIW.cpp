#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

// [[Rcpp::export]]
float cov_loglik (Rcpp::NumericVector u, Rcpp::List params) {
	int N = params["N"];
	int D = params["D"];
	arma::mat S = params["S"];
	
	arma::mat L(D, D, arma::fill::zeros);
	unsigned int k = 0;
	float logDiagL = 1;
	for (int j = 0; j < D; j++) {
		int i = j;
		logDiagL *= u[k];
		for (; i < D; i++) {	 
			L(i, j) = u[k++];
		}
	}
	logDiagL = std::log(logDiagL);
	arma::mat inv_L = inv(L);
	arma::mat t_inv_L = inv_L.t();

	float loglik = -0.5 * N * D * std::log(2*M_PI) - N * logDiagL - 
					0.5 * arma::trace(t_inv_L * inv_L * S);

	return loglik;
}

//TODO make a helper function for lower triag and log diag terms
/*
// [[Rcpp::loglik]]
float test_logprior(Rcpp::NumericVector u, Rcpp::List params) {
	arma::mat Omega = params["Omega"];
	unsigned int nu = params["nu"];
	unsigned int D = params["D"];

	
}
*/

