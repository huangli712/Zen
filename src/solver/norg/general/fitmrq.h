#pragma once

/*
Object for nonlinear least-squares fitting by the Levenberg-Marquardt method, also including
the ability to hold specified parameters at fixed, specified values. Call constructor to bind data
vectors and fitting functions and to input an initial parameter guess. Then call any combination
of hold, free, and fit as often as desired. fit sets the output quantities a, covar, alpha,
and chisq.

Tranlated from Numerical Recipies (in C++) (Third Edition)
by Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China), 2021-02-19
*/

#include "stdpfx.h"
#include "vec.h"
#include "mat.h"
#include "gaussj.h"

template<typename F>
class FitMrq {
public:
	static const Int NDONE = 32;			// default value: 4
	static const Int ITMAX = 2000; //9990;	// default value: 1000
	Int ndat;				// number of data points
	Int ma;					// number of free and fixed fitting parameters to be determined
	Int mfit;				// number of free fitting perameters to be determined
	const VecInt& x;
	const VecReal& y;
	const VecReal& sig;
	Real tol;				// default value: 1.e-3
	// user-supplied function funcs that calculates the nonlinear fitting function and 
	// its derivatives. void (*funcs)(const Int x, const VecReal& a, Real& y, VecReal& dyda);
	const F& funcs;
	VecBool ia;
	VecReal a;				// optimized fitting parameters after fit()
	MatReal covar;			// covariance matrix
	MatReal alpha;
	Real chisq;				// chi_sqr after fit()
private:
	// Used by fit to evaluate the linearized fitting matrix alpha, 
	// and vector beta as in (15.5.8), and to calculate and return chi_sqr.
	Real mrqcof(const VecReal& a, MatReal& alpha, VecReal& beta) const;

	// Expand in storage the covariance matrix covar, so as to take into account 
	// parameters that are being held fixed. (For the latter, return zero covariances.)
	void covsrt(MatReal& covar) const;
public:
	// Constructor. Binds references to the data arrays x_i, y_i, and sig_i, and to a 
	// user-supplied function funks that calculates the nonlinear fitting function and 
	// its derivatives. Also inputs the initial parameters guess a_i (which is copied, 
	// not modified) and an optional convergence tolerance tol_i. Initializes all 
	// parameters as free (not held).
	FitMrq(const VecInt& x_i, const VecReal& y_i, const VecReal& sig_i, const VecReal& a_i,
		const F& funks, const Real tol_i = 1.e-20) :
		ndat(x_i.size()), ma(a_i.size()), x(x_i), y(y_i), sig(sig_i),
		tol(tol_i), funcs(funks), ia(ma, true), a(a_i), covar(ma, ma), alpha(ma, ma) {
	}

	void hold(const Int i, const Real val) { ia[i] = false; a[i] = val; }
	void free(const Int i) { ia[i] = true; }

	// Iterate to reduce the chi_sqr of a fit between a set of data points (x, y) 
	// with individual standard deviations sig, and a nonlinear function that depends on 
	// coefficients a. When chi_sqr is no longer decreasing, set best-fit values for the 
	// parameters a, and chisq = chi_sqr, covar, and alpha. (Parameters held fixed will 
	// return zero covariances.)
	// return (number of actual iterations) if fitting has been implemented successfully
	// return -1 if the number of iterations is exhausted
	// return -2 if fatal error occurs due to singular matrix in gaussj
	Int fit();
};

template<typename F>
Real FitMrq<F>::mrqcof(const VecReal& a, MatReal& alpha, VecReal& beta) const
{
	Int i, j, k, l, m;
	Real ymod, wt, sig2i, dy;
	VecReal dyda(ma);
	for (j = 0; j < mfit; j++) {
		for (k = 0; k <= j; k++) alpha[j][k] = 0.0;
		beta[j] = 0.;
	}
	Real chi_sqr = 0.;
	for (i = 0; i < ndat; i++) {
		funcs(x[i], a, ymod, dyda);
		sig2i = 1.0 / (sig[i] * sig[i]);
		dy = y[i] - ymod;
		for (j = 0, l = 0; l < ma; l++) {
			if (ia[l]) {
				wt = dyda[l] * sig2i;
				for (k = 0, m = 0; m < l + 1; m++)
					if (ia[m]) alpha[j][k++] += wt * dyda[m];	
				beta[j++] += dy * wt;
			}
		}
		chi_sqr += dy * dy * sig2i;
	}
	for (j = 1; j < mfit; j++)
		for (k = 0; k < j; k++) alpha[k][j] = alpha[j][k];
	return chi_sqr;
}

template<typename F>
void FitMrq<F>::covsrt(MatReal& covar) const
{
	Int i, j, k;
	for (i = mfit; i < ma; i++)
		for (j = 0; j < i + 1; j++) covar[i][j] = covar[j][i] = 0.0;
	k = mfit - 1;
	for (j = ma - 1; j >= 0; j--) {
		if (ia[j]) {
			for (i = 0; i < ma; i++) SWAP(covar[i][k], covar[i][j]);
			for (i = 0; i < ma; i++) SWAP(covar[k][i], covar[j][i]);
			k--;
		}
	}
}

template<typename F>
Int FitMrq<F>::fit()
{
	Int j, k, l, iter, done = 0;
	Real alamda = 1.;				// original value: .001
	VecReal atry(ma), beta(ma), da(ma);
	mfit = 0;
	for (j = 0; j < ma; j++) if (ia[j]) mfit++;
	MatReal oneda(mfit, 1), temp(mfit, mfit);
	chisq = mrqcof(a, alpha, beta);
	for (j = 0; j < ma; j++) atry[j] = a[j];
	for (iter = 0; iter < ITMAX; iter++) {
		if (done == NDONE) alamda = 0.;
		for (j = 0; j < mfit; j++) {
			for (k = 0; k < mfit; k++) covar[j][k] = alpha[j][k];
			covar[j][j] = alpha[j][j] * (1.0 + alamda);
			for (k = 0; k < mfit; k++) temp[j][k] = covar[j][k];
			oneda[j][0] = beta[j];
		}
		Int ifgs = gaussj(temp, oneda);
		if (ifgs > 0) {
			return -2;	// ERR("singular matrix in gaussj");
		}
		for (j = 0; j < mfit; j++) {
			for (k = 0; k < mfit; k++) covar[j][k] = temp[j][k];
			da[j] = oneda[j][0];
		}
		if (done == NDONE) {
			covsrt(covar);
			covsrt(alpha);
			return iter;	// fitting has been implemented successfully
		}
		for (j = 0, l = 0; l < ma; l++)
			if (ia[l]) atry[l] = a[l] + da[j++];
		Real chisqtry = mrqcof(atry, covar, da);
		if (ABS(chisq - chisqtry) < MAX(tol, tol * chisq)) done++;
		if (chisqtry < chisq) {
			alamda *= 0.1;
			for (j = 0; j < mfit; j++) {
				for (k = 0; k < mfit; k++) alpha[j][k] = covar[j][k];
				beta[j] = da[j];
			}
			for (l = 0; l < ma; l++) a[l] = atry[l];
			chisq = chisqtry;
		}
		else {
			alamda *= 10;
		}
	}
	return -1;	// ERR("FitMrq too many iterations");
}
