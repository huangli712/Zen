/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China) date 2013 - 2017
	Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2020 - 2022
*/

#ifndef _MILAPACKE_H_
#define _MILAPACKE_H_

#define MKL_Complex8  std::complex<float>
#define MKL_Complex16 std::complex<double>
#include <mkl.h>
#include <mkl_lapacke.h>

// diagonalize Tridiagonal matrix m, eigenvalues are contained in d on output in ascending order
inline Int trd_heevr_qr(Vec<double>& a, Vec<double>& b, Mat<double>& m)
{
	if (!a.size()) return 0;
	lapack_int n = (lapack_int)a.size();
	//double* a = m.p();
	double* d_ = a.p();
	double* e_ = b.p();

	lapack_int info;
	Mat<double> z(n, n);
	Vec<lapack_int> isuppz(2 * n);
	//https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/lapack-routines/lapack-least-squares-and-eigenvalue-problem/lapack-least-square-eigenvalue-problem-computation/symmetric-eigenvalue-problems-lapack-computation/steqr.html
	info = LAPACKE_dsteqr(LAPACK_ROW_MAJOR, 'I', n, d_, e_, z.p(), n);
	m = z;
	return info;
}

// diagonalize Tridiagonal matrix m,
// eigenvalues are contained in d on output in ascending order
// Computes all eigenvalues and, optionally,
// all eigenvectors of a real symmetric tridiagonal matrix using divide and conquer algorithm.
inline Int trd_heevr_vd(Vec<double>& a, Vec<double>& b, Mat<double>& m)
{
	if (!a.size()) return 0;
	lapack_int n = (lapack_int)a.size();
	//double* a = m.p();
	double* d_ = a.p();
	double* e_ = b.p();

	lapack_int info;
	Mat<double> z(n, n);
	Vec<lapack_int> isuppz(2 * n);
	//https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/lapack-routines/lapack-least-squares-and-eigenvalue-problem/lapack-least-squares-eigenvalue-problem-driver/symmetric-eigenvalue-problems-lapack-driver/stevd.html#stevd
	info = LAPACKE_dstevd(LAPACK_ROW_MAJOR, 'V', n, d_, e_, z.p(), n);
	m = z;
	return info;
}

// diagonalize hermitian matrix m
// eigenvalues are contained in d on output in ascending order
// if ifv == 1, the columns of m are replaced by the corresponding eigenvectors on output
inline Int heevr(Mat<double> &m, Vec<double> &d, Int ifv = 1)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(m.nrows(), m.ncols(), m.ncols(), d.size());
#endif
	if (!d.size()) return 0;
	lapack_int n = (lapack_int)d.size();
	double *a = m.p();
	double *w = d.p();
	char jobz = ifv ? 'V' : 'N';
	char range = 'A';
	char uplo = 'U';
	double vl = -1.E+64;
	double vu = 1.E+64;
	double abstol = 0.;
	lapack_int il = 1;
	lapack_int iu = n;
	lapack_int nr;
	lapack_int info;
	Mat<double> z(n, n);
	Vec<lapack_int> isuppz(2 * n);
	info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, jobz, range, uplo, 
		n, a, n, vl, vu, il, iu, abstol, &nr, w, z.p(), n, isuppz.p());
	m = z;
	return info;
}
// diagonalize hermitian matrix m
// eigenvalues are contained in d on output in ascending order
// if ifv == 1, the columns of m are replaced by the corresponding eigenvectors on output
inline Int heevr(Mat<lapack_complex_double> &m, Vec<double> &d, Int ifv = 1)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(m.nrows(), m.ncols(), m.ncols(), d.size());
#endif
	if (!d.size()) return 0;
	lapack_int n = (lapack_int)d.size();
	lapack_complex_double *a = m.p();
	double *w = d.p();
	char jobz = ifv ? 'V' : 'N';
	char range = 'A';
	char uplo = 'U';
	double vl = -1.E+64;
	double vu = 1.E+64;
	double abstol = 0.;
	lapack_int il = 1;
	lapack_int iu = n;
	lapack_int nr;
	lapack_int info;
	Mat<lapack_complex_double> z(n, n);
	Vec<lapack_int> isuppz(2 * n);
	info = LAPACKE_zheevr(LAPACK_ROW_MAJOR, jobz, range, uplo, 
		n, a, n, vl, vu, il, iu, abstol, &nr, w, z.p(), n, isuppz.p());
	m = z;
	return info;
}

// diagonalize hermitian matrix m
// eigenvalues are contained in d on output in ascending order
// if ifv == 1, the columns of m are replaced by the corresponding eigenvectors on output
inline Int heevd(Mat<double> &m, Vec<double> &d, Int ifv = 1)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(m.nrows(), m.ncols(), m.ncols(), d.size());
#endif
	if (!d.size()) return 0;
	lapack_int n = (lapack_int)d.size();
	double *a = m.p();
	double *w = d.p();
	char jobz = ifv ? 'V' : 'N';
	char uplo = 'U';
	return LAPACKE_dsyevd(LAPACK_ROW_MAJOR, jobz, uplo, n, a, n, w);
}
// diagonalize hermitian matrix m
// eigenvalues are contained in d on output in ascending order
// if ifv == 1, the columns of m are replaced by the corresponding eigenvectors on output
inline Int heevd(Mat<lapack_complex_double> &m, Vec<double> &d, Int ifv = 1)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(m.nrows(), m.ncols(), m.ncols(), d.size());
#endif
	if (!d.size()) return 0;
	lapack_int n = (lapack_int)d.size();
	lapack_complex_double *a = m.p();
	double *w = d.p();
	char jobz = ifv ? 'V' : 'N';
	char uplo = 'U';
	return LAPACKE_zheevd(LAPACK_ROW_MAJOR, jobz, uplo, n, a, n, w);
}



// check the singular value decomposition (SVD) results
// A = U �� V^\dagger
// the diagonal elements of �� are the singular values of A; 
// they are real and non-negative, and are returned in descending order.
// s.size() == MIN(a.nrows(), a.ncols()).
// return relative error of SVD: PairReal(erra, erru + errv).
template<typename T>
inline PairReal errsvd(const Mat<T> &a, const Vec<T> &s, const Mat<T> &u, const Mat<T> &v)
{
	lapack_int d = (lapack_int)MIN(a.nrows(), a.ncols());
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ3(a.nrows(), u.nrows(), a.ncols(), v.ncols(), d, s.size());
//	ASSERT_EQ2(u.nrows(), u.ncols(), v.nrows(), v.ncols());
#endif
	Mat<T> vt = v.ct();
	lapack_int m = (lapack_int)a.nrows();
	lapack_int n = (lapack_int)a.ncols();
	Mat<T> u_d(m, d);
	Mat<T> vt_d(d, n);
	Mat<T> Imat = dmat(d, T(1.));
	for_Int (i, 0, m) for_Int (j, 0, d) u_d[i][j] = u[i][j];
	for_Int (i, 0, d) vt_d[i] = vt[i];
	Real erra = RERR(a, u_d * dmatmul(s, vt));
	Real erru = (Imat - u_d.ct() * u_d).norm_avg();
	Real errv = (Imat - vt_d * vt_d.ct()).norm_avg();
	return PairReal(erra, erru + errv);
}

// computes the singular value decomposition of a general rectangular matrix.
// A = U �� V^\dagger
// the diagonal elements of �� are the singular values of A; 
// they are real and non-negative, and are returned in descending order.
// s.size() == MIN(a.nrows(), a.ncols()).
// on output, A is destroyed.
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// this Lapack function may gives wrong results without warning, so we should always check the results.       !!!!!!!!!!!!!!
// it give a wrong result for an A of dimension 1024 * 512 with random elements in [-0.5, 0.5]                !!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
inline Int gesvd(Mat<double> &a, Vec<double> &s, Mat<double> &u, Mat<double> &v)
{
	lapack_int mind = (lapack_int)MIN(a.nrows(), a.ncols());
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ3(a.nrows(), u.nrows(), a.ncols(), v.ncols(), mind, s.size());
	ASSERT_EQ2(u.nrows(), u.ncols(), v.nrows(), v.ncols());
#endif
	if (!a.size()) return 0;

	char jobu = 'A';
	char jobvt = 'A';
	lapack_int m = (lapack_int)a.nrows();
	lapack_int n = (lapack_int)a.ncols();
	lapack_int lda = n;
	lapack_int ldu = m;
	lapack_int ldvt = n;
	Mat<double> vt(n, n);
	Vec<double> superb(mind);
	Mat<double> atmp = a;

	lapack_int info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, jobu, jobvt, m, n, a.p(), lda, s.p(), u.p(), ldu, vt.p(), ldvt, superb.p());
	v = vt.ct();

	PairReal err = errsvd(atmp, s, u, v);
	double erra = err.first;
	double erru = err.second;
	double eps = eps_acc_def(MAX(a.nrows(), a.ncols()));
	double epsbig = eps_acc(MAX(a.nrows(), a.ncols()), SQRT(eps_Real));
	if (info || erra > epsbig || erru > epsbig) ERR(NAV6(m, n, info, erra, erru, epsbig));
	if (erra > eps || erru > eps) WRN(NAV6(m, n, info, erra, erru, eps));

	return info;
}
// computes the singular value decomposition of a general rectangular matrix.
// A = U �� V^\dagger
// the diagonal elements of �� are the singular values of A; 
// they are real and non-negative, and are returned in descending order.
// s.size() == MIN(a.nrows(), a.ncols()).
// on output, A is destroyed.
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// this Lapack function may gives wrong results without warning, so we should always check the results.       !!!!!!!!!!!!!!
// it give a wrong result for an A of dimension 1024 * 512 with random elements in [-0.5, 0.5]                !!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
inline Int gesvd(Mat<double> &a, Vec<double> &s)
{
	lapack_int mind = (lapack_int)MIN(a.nrows(), a.ncols());
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(mind, s.size());
#endif
	if (!a.size()) return 0;

	char jobu = 'N';
	char jobvt = 'N';
	lapack_int m = (lapack_int)a.nrows();
	lapack_int n = (lapack_int)a.ncols();
	lapack_int lda = n;
	lapack_int ldu = m;
	lapack_int ldvt = n;
	Mat<double> u;
	Mat<double> vt;
	Vec<double> superb(mind);
	lapack_int info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, jobu, jobvt, m, n, a.p(), lda, s.p(), u.p(), ldu, vt.p(), ldvt, superb.p());
	if (info) ERR(NAV3(m, n, info));
	return info;
}


// computes the singular value decomposition of a general rectangular matrix using a divide and conquer method.
// A = U �� V^\dagger
// the diagonal elements of �� are the singular values of A; 
// they are real and non-negative, and are returned in descending order.
// s.size() == MIN(a.nrows(), a.ncols()).
// on output, A is destroyed.
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// this Lapack function may gives wrong results without warning, so we should always check the results.       !!!!!!!!!!!!!!
// it give a wrong result for an A of dimension 1024 * 512 with random elements in [-0.5, 0.5]                !!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
inline Int gesdd(Mat<double> &a, Vec<double> &s, Mat<double> &u, Mat<double> &v)
{
	lapack_int mind = (lapack_int)MIN(a.nrows(), a.ncols());
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ3(a.nrows(), u.nrows(), a.ncols(), v.ncols(), mind, s.size());
	ASSERT_EQ2(u.nrows(), u.ncols(), v.nrows(), v.ncols());
#endif
	if (!a.size()) return 0;
	
	char jobz = 'A';
	lapack_int m = (lapack_int)a.nrows();
	lapack_int n = (lapack_int)a.ncols();
	lapack_int lda = n;
	lapack_int ldu = m;
	lapack_int ldvt = n;
	Mat<double> vt(n, n);
	Mat<double> atmp = a;

	lapack_int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, jobz, m, n, a.p(), lda, s.p(), u.p(), ldu, vt.p(), ldvt);
	v = vt.ct();

	PairReal err = errsvd(atmp, s, u, v);
	double erra = err.first;
	double erru = err.second;
	double eps = eps_acc_def(MAX(a.nrows(), a.ncols()));
	double epsbig = eps_acc(MAX(a.nrows(), a.ncols()), SQRT(eps_Real));
	if (info || erra > epsbig || erru > epsbig) ERR(NAV6(m, n, info, erra, erru, epsbig));
	if (erra > eps || erru > eps) WRN(NAV6(m, n, info, erra, erru, eps));

	return info;
}
// computes the singular value decomposition of a general rectangular matrix using a divide and conquer method.
// A = U �� V^\dagger
// the diagonal elements of �� are the singular values of A; 
// they are real and non-negative, and are returned in descending order.
// s.size() == MIN(a.nrows(), a.ncols()).
// on output, A is destroyed.
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// this Lapack function may gives wrong results without warning, so we should always check the results.       !!!!!!!!!!!!!!
// it give a wrong result for an A of dimension 1024 * 512 with random elements in [-0.5, 0.5]                !!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
inline Int gesdd(Mat<double> &a, Vec<double> &s)
{
	lapack_int mind = (lapack_int)MIN(a.nrows(), a.ncols());
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(mind, s.size());
#endif
	if (!a.size()) return 0;

	char jobz = 'N';
	lapack_int m = (lapack_int)a.nrows();
	lapack_int n = (lapack_int)a.ncols();
	lapack_int lda = n;
	lapack_int ldu = m;
	lapack_int ldvt = n;
	Mat<double> u;
	Mat<double> vt;
	lapack_int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, jobz, m, n, a.p(), lda, s.p(), u.p(), ldu, vt.p(), ldvt);
	if (info) ERR(NAV3(m, n, info));
	return info;
}

// computes the singular value decomposition using a preconditioned Jacobi SVD method.
// A = U �� V^\dagger
// A.nrows() >= A.ncols()
// on output, A is NOT destroyed.
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// this Lapack function is about three times slower than other SVD functions in Lapack                        !!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// their are bugs in this Lapack function when m > n (e.g., m = 4, n = 1) and jobu = 'F'                      !!!!!!!!!!!!!!
// it even gives a wrong info and warning, that is "info = -18, but it says parameter 17 is incorrect.        !!!!!!!!!!!!!!
// it turns out that it gives correct results when we set jobu = 'U'.                                         !!!!!!!!!!!!!!
// so the function can be used but it is dangerous. We should always check the results.                       !!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
inline Int gejsv(const Mat<double> &a, Vec<double> &sva, Mat<double> &u, Mat<double> &v)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ3(a.nrows(), u.nrows(), a.ncols(), v.ncols(), v.ncols(), sva.size());
	ASSERT_EQ2(u.nrows(), u.ncols(), v.nrows(), v.ncols());
	ASSERT_GE(a.nrows(), a.ncols());
#endif
	if (!a.size()) return 0;

	char joba = 'F';
	char jobu = 'U';
	char jobv = 'V';
	char jobr = 'R';
	char jobt = 'T';
	char jobp = 'P';
	lapack_int m = (lapack_int)a.nrows();
	lapack_int n = (lapack_int)a.ncols();
	lapack_int lda = n;
	lapack_int ldu = m;
	lapack_int ldv = n;
	Vec<double> stat(7);
	Vec<lapack_int> istat(3);
	Mat<double> atmp = a;
	
	lapack_int info = LAPACKE_dgejsv(LAPACK_ROW_MAJOR, joba, jobu, jobv, jobr, jobt, jobp, 
		m, n, a.p(), lda, sva.p(), u.p(), ldu, v.p(), ldv, stat.p(), istat.p());

	PairReal err = errsvd(atmp, sva, u, v);
	double erra = err.first;
	double erru = err.second;
	double eps = eps_acc_def(MAX(a.nrows(), a.ncols()));
	double epsbig = eps_acc(MAX(a.nrows(), a.ncols()), SQRT(eps_Real));
	if (info || erra > epsbig || erru > epsbig) ERR(NAV6(m, n, info, erra, erru, epsbig));
	if (erra > eps || erru > eps) WRN(NAV6(m, n, info, erra, erru, eps));

	return info;
}
// computes the singular value decomposition using a preconditioned Jacobi SVD method.
// A = U �� V^\dagger
// A.nrows() >= A.ncols()
// on output, A is NOT destroyed.
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// this Lapack function may gives wrong results without warning, so we should always check the results.       !!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
inline Int gejsv(const Mat<double> &a, Vec<double> &sva)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(a.ncols(), sva.size());
	ASSERT_GE(a.nrows(), a.ncols());
#endif
	if (!a.size()) return 0;

	char joba = 'F';
	char jobu = 'N';
	char jobv = 'N';
	char jobr = 'R';
	char jobt = 'T';
	char jobp = 'P';
	lapack_int m = (lapack_int)a.nrows();
	lapack_int n = (lapack_int)a.ncols();
	lapack_int lda = n;
	lapack_int ldu = m;
	lapack_int ldv = n;
	Vec<double> stat(7);
	Vec<lapack_int> istat(3);
	Mat<double> u(m, m);
	Mat<double> v(n, n);
	
	lapack_int info = LAPACKE_dgejsv(LAPACK_ROW_MAJOR, joba, jobu, jobv, jobr, jobt, jobp, 
		m, n, a.p(), lda, sva.p(), u.p(), ldu, v.p(), ldv, stat.p(), istat.p());
	if (info) ERR(NAV3(m, n, info));
	return info;
}

// computes the singular value decomposition of a general rectangular matrix.
// A = U �� V^\dagger
// the diagonal elements of �� are the singular values of A; 
// they are real and non-negative, and are returned in descending order.
// s.size() == MIN(a.nrows(), a.ncols()).
// on output, A is destroyed.
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// this Lapack function may gives wrong results without warning, so we should always check the results.       !!!!!!!!!!!!!!
// it give a wrong result for an A of dimension 1024 * 512 with random elements in [-0.5, 0.5]                !!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
inline Int gesvd(MatCmplx &a, Vec<double> &s, MatCmplx &u, MatCmplx &v)
{
	lapack_int mind = (lapack_int)MIN(a.nrows(), a.ncols());
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ3(a.nrows(), u.nrows(), a.ncols(), v.ncols(), mind, s.size());
	ASSERT_EQ2(u.nrows(), u.ncols(), v.nrows(), v.ncols());
#endif
	if (!a.size()) return 0;

	char jobu = 'A';
	char jobvt = 'A';
	lapack_int m = (lapack_int)a.nrows();
	lapack_int n = (lapack_int)a.ncols();
	lapack_int lda = n;
	lapack_int ldu = m;
	lapack_int ldvt = n;
	MatCmplx vt(n, n);
	Vec<double> superb(mind);
	MatCmplx atmp = a;

	lapack_int info = LAPACKE_zgesvd(LAPACK_ROW_MAJOR, jobu, jobvt, m, n, a.p(), lda, s.p(), u.p(), ldu, vt.p(), ldvt, superb.p());
	v = vt.ct();

	/*PairReal err = errsvd(atmp, s, u, v);
	double erra = err.first;
	double erru = err.second;
	double eps = eps_acc_def(MAX(a.nrows(), a.ncols()));
	double epsbig = eps_acc(MAX(a.nrows(), a.ncols()), SQRT(eps_Real));
	if (info || erra > epsbig || erru > epsbig) ERR(NAV6(m, n, info, erra, erru, epsbig));
	if (erra > eps || erru > eps) WRN(NAV6(m, n, info, erra, erru, eps));*/

	return info;
}

// calculate inverse of matrix ain
inline Mat<double> matinvlu(const Mat<double> &ain)
{
	Mat<double> a = ain;
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(a.nrows(), a.ncols());
#endif
	lapack_int n = (lapack_int)a.nrows();
	lapack_int lda = n;
	Vec<lapack_int> ipiv(n);
	lapack_int info;

	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a.p(), lda, ipiv.p());
#ifdef _ASSERTION_
	ASSERT_EQ(info, 0);
#endif
	info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, a.p(), lda, ipiv.p());
#ifdef _ASSERTION_
	ASSERT_EQ(info, 0);
#endif
	return a;
}
// calculate inverse of matrix ain
inline Mat<lapack_complex_double> matinvlu(const Mat<lapack_complex_double> &ain)
{
	Mat<lapack_complex_double> a = ain;
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(a.nrows(), a.ncols());
#endif
	lapack_int n = (lapack_int)a.nrows();
	lapack_int lda = n;
	Vec<lapack_int> ipiv(n);
	lapack_int info;

	info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, a.p(), lda, ipiv.p());
#ifdef _ASSERTION_
	ASSERT_EQ(info, 0);
#endif
	info = LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, a.p(), lda, ipiv.p());
#ifdef _ASSERTION_
	ASSERT_EQ(info, 0);
#endif
	return a;
}

#endif /* _MILAPACKE_H_ */
