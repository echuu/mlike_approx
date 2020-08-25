#include "fast_psi.h"

using namespace Rcpp;

// [[Rcpp::export]]
float psi (NumericVector& u, List& prior) {
	NumericVector y = prior["y"];
	NumericMatrix X = prior["X"];
	const unsigned int N = prior["N"];
	const unsigned int D = prior["D"]; 
	const float sigmasq = prior["sigmasq"];	
	const float tau = prior["tau"];	
	

	unsigned int nrow = X.nrow();
	unsigned int len = u.size(); 
	float sum = 0;
	for (unsigned int i = 0; i < nrow; i++) {
		float accum = 0;
		for (unsigned int j = 0; j < len; j++) {
			accum += u[j] * X(i, j);
		}
		sum += pow(y[i] - accum, 2);
	}


	float loglik = -0.5 * N * std::log(2 * M_PI * sigmasq) - 
					1.0 / (2.0 * sigmasq) * sum; 
	
	float logTN = 0;
	const float term1 = std::log(1 / std::sqrt((sigmasq  * 2.0 * M_PI) / tau));
	const float term2 = 2 * (sigmasq / tau);
	for (unsigned int i = 0; i < len; i++) {
		logTN += dnorm_log(u[i], 0, term2, term1);
	}


	logTN += (D * std::log(2));
	 return -loglik - logTN;

}

