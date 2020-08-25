#ifndef FAST_PSI
#define FAST_PSI
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Rcpp.h>
#include <cmath>

/*
Calculates normal distribution log 

Args:
	x: 
	mean: mean
	term1: std::log((1.0/(std::sqrt(2.0 * M_PI * sigmasq))));	this is a constant term to prevent recomputation
	sigma2: 2.0 * sigmasq
*/
inline float dnorm_log (const float& x, const float& mean, const float& sigma2, const float& term1) {
	return term1 - std::pow(x-mean, 2) / sigma2; 
}
#endif

