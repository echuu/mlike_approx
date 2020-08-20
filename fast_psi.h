#ifndef FAST_PSI
#define FAST_PSI
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Rcpp.h>
#include <cmath>

float psi(Rcpp::NumericVector u, Rcpp::List prior);

#endif

