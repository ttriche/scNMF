#define ARMA_NO_DEBUG

#ifdef _OPENMP
#include <omp.h>
#endif

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <Rcpp.h>
#include <R.h>

//using namespace Rcpp;
using namespace arma;

/* ****************************************

This code is heavily modified from the NNLM::nnmf GitHub package
An emphasis was placed on increasing computation speed for specific use cases
The method implements the equivalent of:

NNLM::nnmf(...method = "scd", loss = "mse", alpha = c(0,0,0), beta = c(0,0,0), mask = NULL, init = NULL)

All irrelevant conditional statements (given above criteria) were removed
MKL error calculations were removed
pragma omp critical section was removed from for loop in update function
MSE loss is calculated in a parallel manner
Additional speed gains may be obtained by sparse matrix support
However, when zeros aren't missing values, sparse operations will not offer dramatic speed gains

*************************************** */

Rcpp::List c_nnmf(const arma::mat &A,
				  const unsigned int k,
				  const unsigned int max_iter,
				  const double rel_tol,
				  const int n_threads,
				  const int verbose,
				  const unsigned int inner_max_iter,
				  const double inner_rel_tol,
				  const int trace);

Rcpp::List c_project(const arma::mat &A,
					 const mat &W,
					 const int n_threads,
					 const unsigned int inner_max_iter,
					 const double inner_rel_tol);

// calculates least squares with sequential coordinate descent
vec scd_ls_update(subview_col<double> Hj,
				  const arma::mat &WtW, vec &mu,
				  const unsigned int &inner_max_iter,
				  const double &inner_rel_tol)
{
	double tmp;
	double etmp = 0;
	double inner_rel_err = 1 + inner_rel_tol;
	for (unsigned int t = 0; t < inner_max_iter && inner_rel_err > inner_rel_tol; t++)
	{
		inner_rel_err = 0;
		for (unsigned int k = 0; k < WtW.n_cols; k++)
		{
			tmp = Hj(k) - mu(k) / WtW(k, k);
			if (tmp < 0)
				tmp = 0;
			if (tmp != Hj(k))
				mu += (tmp - Hj(k)) * WtW.col(k);
			else
				continue;
			etmp = 2 * std::abs(Hj(k) - tmp) / (tmp + Hj(k) + 1e-16);
			if (etmp > inner_rel_err)
				inner_rel_err = etmp;
			Hj(k) = tmp;
		}
	}
	return (Hj);
}

// updates H given X and A (updates W given Ht and At)
arma::mat ls_update(arma::mat &H,
					const arma::mat &Wt,
					const arma::mat &A,
					unsigned int inner_max_iter,
					double inner_rel_tol,
					int n_threads)
{
	arma::mat WtW = Wt * Wt.t();
	WtW.diag() += 1e-16;

// multithreaded computation of least squares for each column in A with sequential coordinate descent
#pragma omp parallel for num_threads(n_threads) schedule(dynamic)
	for (unsigned int j = 0; j < A.n_cols; j++)
	{
		vec mu = WtW * H.col(j) - Wt * A.col(j);
		H.col(j) = scd_ls_update(H.col(j), WtW, mu, inner_max_iter, inner_rel_tol);
	}
	return (H);
}

// calculates mean squared error given A, W, H, and a number of threads for parallelization
double calc_mse(const arma::mat &A,
				const arma::mat &W,
				const arma::mat &H,
				int n_threads)
{
	double mse_err = 0;
	arma::mat WtH = W.t() * H;
	unsigned int A_ncols = A.n_cols;

#pragma omp parallel for num_threads(n_threads) schedule(dynamic)
	for (unsigned int c = 0; c < A_ncols; c++)
	{
		mse_err += mean(square(A.col(c) - WtH.col(c)));
	}
	return (mse_err / A_ncols);
}

//[[Rcpp::export]]
Rcpp::List c_nnmf(const arma::mat &A,
				  const unsigned int k,
				  const unsigned int max_iter,
				  const double rel_tol,
				  const int n_threads,
				  const int verbose,
				  const unsigned int inner_max_iter,
				  const double inner_rel_tol,
				  const int trace)
{

	unsigned int err_len = (unsigned int)std::ceil(double(max_iter)) + 1;
	double rel_err = rel_tol + 1;
	vec mse_err(err_len);

	if (verbose > 0)
		Rprintf("\n%10s | %10s | %10s\n", "Iteration", "MSE", "Rel. Err.\n-----------------------------------------------");

	// randomly initialize W and H
	arma::mat W;
	arma::mat H;
	W = W.randu(k, A.n_rows);
	H = H.randu(k, A.n_cols);
	W *= 0.01;
	H *= 0.01;

	// loop for alternating least squares updating of W and H with error checking
	for (unsigned int i = 0; i < max_iter && std::abs(rel_err) > rel_tol; i++)
	{
		Rcpp::checkUserInterrupt();

		// update W and then H
		W = ls_update(W, H, A.t(), inner_max_iter, inner_rel_tol, n_threads);
		H = ls_update(H, W, A, inner_max_iter, inner_rel_tol, n_threads);

		// calculate error if trace conditions are met
		if (i % trace == 0 || ((i + 1) % trace == 0))
			mse_err(i) = calc_mse(A, W, H, n_threads);

		// update rel_err based on err of iteration i-1
		if (i % trace == 0 && i > 3)
		{
			rel_err = (mse_err(i - 1) - mse_err(i)) / (mse_err(i - 1) + mse_err(i) + 1e-16);
			if (verbose > 0)
				Rprintf("%10d | %10.4f | %10.g\n", i + 1, mse_err(i), rel_err);
		}
	}

	if (rel_err > rel_tol)
		Rcpp::warning("Target tolerance not reached within max.iter iterations. Try a larger max.iter.");

	return Rcpp::List::create(Rcpp::Named("W") = W.t(), Rcpp::Named("H") = H);
}

//[[Rcpp::export]]
Rcpp::List c_project(const arma::mat &A,
					 const arma::mat &W,
					 const int n_threads,
					 const unsigned int inner_max_iter,
					 const double inner_rel_tol)
{
	// randomly initialize H
	arma::mat H;
	H = H.randu(W.n_cols, A.n_cols);
	H *= 0.01;

	// given A and W, calculate H. This is a single instance of the minimization problem (A - WH)^2
	H = ls_update(H, W.t(), A, inner_max_iter, inner_rel_tol, n_threads);

	return Rcpp::List::create(Rcpp::Named("H") = H);
}
