/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2015-08-28
*/

#ifndef _LUDCMP_H_
#define _LUDCMP_H_

#include "stdpfx.h"
#include "vec.h"
#include "mat.h"

template<typename T>
// P * A = L * U, or A = P^-1 * L * U, 
// where P is a rowwise permutation, L is lower triangular (has elements only on the diagonal and below) and 
// U is upper triangular (has elements only on the diagonal and above).
// the diagonal unity elements L_{ii} are not stored at all.
// pm is an output vector that records the row permutation effected by the partial pivoting.
// d is output as +1/-1 depending on whether the number of row interchanges was even or odd, respectively.
// pivot technique is used if pv == 1, otherwise not and P == 1 is imposed.
// L and U are stored in matrix lu.
class LUdcmp {
private:
	const Mat<T> &a;		// a references A, used only by improve
	Int n;					// n = a.nrows() = a.ncols()
	Int pv;					// if use pivoting
	Int sc;					// singularity code, 0 for non-singular, 1 for singular with zero-rows, 2 for singular after nontrivial examination
	Mat<T> lu;				// stores L and U
	Real d;					// parity, used by det()
	VecInt pm;				// stores the permutation
public:
	inline LUdcmp(const Mat<T> &a_i, Int pivot = 1);
	inline void solve(const Vec<T> &b, Vec<T> &x) const;
	inline void solve(const Mat<T> &b, Mat<T> &x) const;
	inline void improve(const Vec<T> &b, Vec<T> &x) const;
	inline T det() const;
	inline void inverse(Mat<T> &ainv) const;
	inline Int singular() const { return sc; }
	inline Real parity() const { return d; }
	inline Int dim() const { return n; }
	inline const Vec<T> &operator[](Idx i) const { return lu[i]; }
	template<typename Tin> inline void permute(Mat<Tin> &a) const;
	template<typename Tin> inline void permute(Vec<Tin> &v) const;
};

// related functions

template<typename T> T determinant(const Mat<T> &a) { return LUdcmp<T>(a).det(); }

// function definitions

template<typename T>
// this constructor has done a LU composition for matrix a, L and U are stored in matrix lu
LUdcmp<T>::LUdcmp(const Mat<T> &a_i, Int pivot): a(a_i), n(a.nrows()), pv(pivot), sc(0), lu(a), pm(n)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(a.nrows(), a.ncols());
#endif
	const Real TINY=1.E-80;
	Int i,imax,j,k;
	Real big,mgn;
	Vec<T> vtmp(n);
	VecReal vv(n);
	for (i=0;i<n;i++) {
		big=0.;
		for (j=0;j<n;j++)
			if ((mgn = ABS(lu[i][j])) > big) big = mgn;
		if (big == 0.) {
			lu[i][i] = big = TINY;
			sc = 1;						// singular matrix in LUdcmp, a zero row in the target matrix.
		}
		vv[i]=1./big;
	}
	d=1.;
	for (k=0;k<n;k++) {
		imax=k;
		big=vv[k]*ABS(lu[k][k]);
		for (i=k+1;i<n;i++) {
			mgn=vv[i]*ABS(lu[i][k]);
			if (pv == 1 && mgn > big) {
				big=mgn;
				imax=i;
			}
		}
		if (k != imax) {
			vtmp = lu[imax];
			lu[imax] = lu[k];
			lu[k] = vtmp;
			d = -d;
			vv[imax]=vv[k];
		}
		pm[k]=imax;
		if (lu[k][k] == 0.) {
			lu[k][k]=TINY;
			sc = 2;				// singular matrix in LUdcmp, nontrivial examination.
		}
		for (i=k+1;i<n;i++) {
			T tmp=lu[i][k] /= lu[k][k];
			for (j=k+1;j<n;j++) lu[i][j] -= tmp*lu[k][j];
		}
	}
}

template<typename T>
// solve A x = b
// on input b and x may reference the same vector, 
// in which case the solution overwrites the input.
void LUdcmp<T>::solve(const Vec<T> &b, Vec<T> &x) const
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(b.size(), n, x.size(), n);
#endif
	Int i,ii=0,ip,j;
	T sum;
	x = b;
	for (i=0;i<n;i++) {
		ip=pm[i];
		sum=x[ip];
		x[ip]=x[i];
		if (ii != 0)
			for (j=ii-1;j<i;j++) sum -= lu[i][j]*x[j];
		else if (sum != 0.)
			ii=i+1;
		x[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=x[i];
		for (j=i+1;j<n;j++) sum -= lu[i][j]*x[j];
		x[i]=sum/lu[i][i];
	}
}

template<typename T>
// solve A X = B
// on input b and x may reference the same matrix, 
// in which case the solution overwrites the input.
void LUdcmp<T>::solve(const Mat<T> &b, Mat<T> &x) const
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ3(b.nrows(), n, x.nrows(), n, b.ncols(), x.ncols());
#endif
	int i,j,m=b.ncols();
	Vec<T> xx(n);
	for (j=0;j<m;j++) {
		for (i=0;i<n;i++) xx[i] = b[i][j];
		solve(xx,xx);
		for (i=0;i<n;i++) x[i][j] = xx[i];
	}
}

template<typename T>
// improve a solution vector x of A x = b, i.e., x will be more accurate after this function call.
// you can call this function several times in succession if you want. 
// Unless you are starting quite far from the true solution, one call is generally enough; 
// but a second call to verify convergence can be reassuring.
void LUdcmp<T>::improve(const Vec<T> &b, Vec<T> &x) const
{
	Vec<T> r(n);
	for_Int (i, 0, n) {
		LReal sdp = -b[i];
		for_Int (j, 0, n) sdp += LReal(a[i][j]) * LReal(x[j]);
		r[i] = sdp;
	}
	solve(r, r);
	x -= r;
}

template<typename T>
T LUdcmp<T>::det() const
{
	const Real TINY=1.E-80;
	T dd = d;
	for (Int i=0;i<n;i++) dd *= lu[i][i];
	if (sc) {
		if (ABS(dd) > SQRT(TINY)) WRN(NAV3(sc, dd, TINY))
		else dd = 0.;
	}
	return dd;
}

template<typename T>
// return the inverse of matrix A in ainv
void LUdcmp<T>::inverse(Mat<T> &ainv) const
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(ainv.nrows(), n, ainv.ncols(), n);
#endif
	ainv = dmat(n, 1.);
	solve(ainv, ainv);
}

template<typename T>
template<typename Tin>
// rowwise permutation of a, return P * a
void LUdcmp<T>::permute(Mat<Tin> &a) const
{
	Vec<T> tmp(n);
	for_Int (i, 0, n) {
		if (i != pm[i]) {
			tmp = a[i];
			a[i] = a[pm[i]];
			a[pm[i]] = tmp;
		}
	}
}

template<typename T>
template<typename Tin>
// rowwise permutation of v, return P * v, where P * A = L * U
void LUdcmp<T>::permute(Vec<Tin> &v) const
{
	T tmp;
	for_Int (i, 0, n) {
		if (i != pm[i]) {
			tmp = v[i];
			v[i] = v[pm[i]];
			v[pm[i]] = tmp;
		}
	}
}


#endif /* _LUDCMP_H_ */
