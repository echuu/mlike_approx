#include <RcppArmadillo.h>
#include <cmath>


typedef unsigned int u_int;


// [[Rcpp::depends(RcppArmadillo)]]
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
float create_loTriag(arma::mat& L, Rcpp::NumericVector& u);
inline float logmultigamma (unsigned int p, float a);
float cov_loglik (Rcpp::NumericVector& u, Rcpp::List& params); 
float cov_logprior(Rcpp::NumericVector& u, Rcpp::List& params); 


// [[Rcpp::export]]
inline float psi(Rcpp::NumericVector& u, Rcpp::List& params) {
	return -cov_logprior(u, params); 
}

float cov_loglik (Rcpp::NumericVector& u, Rcpp::List& params) {
	int N = params["N"];
	int D = params["D"];
	arma::mat S = params["S"];
	
	arma::mat L(D, D, arma::fill::zeros);

	float logDiagL = create_loTriag(L, u);
	arma::mat inv_L = inv(L);
	arma::mat t_inv_L = inv_L.t();

	float loglik = -0.5 * N * D * std::log(2*M_PI) - N * logDiagL - 
					0.5 * arma::trace(t_inv_L * inv_L * S);

	return loglik;
}

// [[Rcpp::export]]
arma::mat vec2chol(Rcpp::NumericVector& u, u_int D) {

	arma::mat Lt(D, D, arma::fill::zeros);
	u_int p = 0;
	for (u_int c = 0; c < D; c++) {
		for (u_int r = 0; r <= c; r++) {
			Lt(r, c) = u[p++];	
		}
	}
	return Lt;
}



// [[Rcpp::export]]
Rcpp::NumericVector chol2vec(arma::mat& M, u_int D) {

	u_int k = 0; 
	u_int D_0 = D * (D + 1) / 2;
	Rcpp::NumericVector u(D_0);
	for (u_int c = 0; c < D; c++) {
		for (u_int r = 0; r <= c; r++) {
			u[k++] = M(r, c);	
		}
	}
	return u;
}



// [[Rcpp::export]]
Rcpp::NumericVector grad(Rcpp::NumericVector& u, Rcpp::List& params) {
	arma::vec nu     = params["nu"];
	arma::vec xi     = params["xi"];
	arma::mat V      = params["V"];
	u_int D   = params["D"];
	u_int D_0 = params["D_0"];

	arma::mat Lt = vec2chol(u, D);

	arma::mat diag_terms(D, D, arma::fill::zeros);
	for (u_int i = 0; i < D; i++) { diag_terms(i,i) = xi[i] / Lt(i,i); }
	
	arma::mat grad(D, D, arma::fill::zeros);
	grad = - diag_terms + Lt * V;

	Rcpp::NumericVector grad_vec = chol2vec(grad, D);

	return grad_vec;
}

// [[Rcpp::export]]
arma::mat hess(Rcpp::NumericVector& u, Rcpp::List& params) {

	arma::vec nu      = params["nu"];
	arma::vec xi      = params["xi"];
	arma::mat V       = params["V"];
	unsigned int D    = params["D"];
	unsigned int D_0  = params["D_0"];
	arma::mat ind_mat = params["t_ind"];

	// initialize the output matrix
	arma::mat H(D_0, D_0, arma::fill::zeros);

	// reconstruct upper cholesky factor -- to be moved into another function
	arma::mat Lt = vec2chol(u, D);

	unsigned int i, j, k, l, r, c;
	arma::mat test(D_0, 2, arma::fill::zeros);
	for (r = 0; r < D_0; r++) {
		i = ind_mat(r, 0) - 1; // subtract one to account for 0-index
		j = ind_mat(r, 1) - 1;
		c = r;
		while (c < D_0) {

			k = ind_mat(c, 0) - 1; // row of 2nd order partial
			l = ind_mat(c, 1) - 1; // col of 2nd order partial

			if (i != k) {
				H(r,c) = 0;
			} else if (i == j && j == k && k == l) {
				H(r,c) = xi(i) / std::pow(Lt(i,i), 2) + V(i,i);
			} else if (i != j && k == i && l >= j) {
				H(r,c) = V(l,j);
			}
			c++;
		} // end inner while loop

	} // end for() populating hessian matrix

	H = 0.5 * (H + H.t());

	return H;
}



float cov_logprior(Rcpp::NumericVector& u, Rcpp::List& params) {

	arma::mat V     = params["V"];
	unsigned int b  = params["b"];
	unsigned int D  = params["D"];
	arma::vec nu    = params["nu"];
	arma::mat Lt(D, D, arma::fill::zeros);

	unsigned int k = 0;
	for (unsigned int c = 0; c < D; c++) {
		for (unsigned int r = 0; r <= c; r++) {
			Lt(r, c) = u[k++];	
		}
	}

	float logprior = 0;
	for (unsigned int i = 0; i < D; i++) {
		logprior += (b + nu(i) - 1) * std::log(Lt(i,i));
	}
	logprior += - 0.5 * arma::trace(Lt.t() * Lt * V) + D * std::log(2);

	return logprior;
}

/* Populars lower diagonal of N x N matrix L with values from u,
	params:
		L: N x N arma matrix
		u: numeric vector
	returns:
		sum of log diag values
*/
float create_loTriag(arma::mat& L, Rcpp::NumericVector& u) {
	if (L.n_rows != L.n_cols)
		throw std::invalid_argument("Unequal dimensions");

	unsigned int D = L.n_rows;
	unsigned int k = 0;	
	float logDiagL = 1;

	for (unsigned int j = 0; j < D; j++) {
		logDiagL *= u[k];
		for (unsigned i = j; i < D; i++) {
			L(i, j) = u[k++];	
		}
	}
	return std::log(logDiagL);
}

/* log multivariate gamma function
	params:
*/
inline float logmultigamma (unsigned int p, float a) {
	float f = 0.25 * p * (p - 1) * std::log(M_PI);
	for (unsigned int i = 1; i <= p; i++)
		f += lgamma(a+0.5-(0.5 * i));
	return f;
}

