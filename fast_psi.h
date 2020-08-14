#ifndef FAST_PSI
#define FAST_PSI
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>

#define MAT_TYPE Eigen::MatrixXd

float dnorm_log(float x, float mean, float sd);

float psi(Rcpp::NumericVector u, Rcpp::List prior);

#endif

