/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2013 - 2017
*/

#ifndef _MILAPACK_H_
#define _MILAPACK_H_

#define MKL_Complex8  std::complex<float>
#define MKL_Complex16 std::complex<double>
#include <mkl.h>

// diagonalize matrix m, eigenvalues are contained in d on output in ascending order
// if ifv == 1, eigenvectors will be found and contained in m on output
inline Int heevr(Mat<double> &m, Vec<double> &d, Int ifv = 1)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(m.nrows(), m.ncols(), m.ncols(), d.size());
#endif
	if (!d.size()) return 0;
	MKL_INT n = d.size();
	double *a = m.p();
	double *w = d.p();
	char jobz = ifv ? 'V' : 'N';
	char range = 'A';
	char uplo = 'U';
	double vl = -1.E+64;
	double vu = 1.E+64;
	double abstol = 0.;
	MKL_INT il = 1;
	MKL_INT iu = n;
	MKL_INT nr;
	MKL_INT info;
	Mat<double> z(n, n);
	Vec<MKL_INT> isuppz(2 * MAX(1, n));
	MKL_INT lwork = MAX(1, 26 * n);
	Vec<double> work(lwork);
	MKL_INT liwork = MAX(1, 10 * n);
	Vec<MKL_INT> iwork(liwork);
	dsyevr(&jobz, &range, &uplo, &n, a, &n, &vl, &vu, &il, &iu, &abstol, &nr, w, z.p(), &n, isuppz.p(), 
		work.p(), &lwork, iwork.p(), &liwork, &info);
	m = z.tr();
	return info;
}

// diagonalize matrix m, eigenvalues are contained in d on output in ascending order
// if ifv == 1, eigenvectors will be found and contained in m on output
inline Int heevr(Mat<MKL_Complex16> &m, Vec<double> &d, Int ifv = 1)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(m.nrows(), m.ncols(), m.ncols(), d.size());
#endif
	if (!d.size()) return 0;
	MKL_INT n = d.size();
	MKL_Complex16 *a = m.p();
	double *w = d.p();
	char jobz = ifv ? 'V' : 'N';
	char range = 'A';
	char uplo = 'U';
	double vl = -1.E+64;
	double vu = 1.E+64;
	double abstol = 0.;
	MKL_INT il = 1;
	MKL_INT iu = n;
	MKL_INT nr;
	MKL_INT info;
	Mat<MKL_Complex16> z(n, n);
	Vec<MKL_INT> isuppz(2 * MAX(1, n));
	MKL_INT lwork = MAX(1, 2 * n);
	Vec<MKL_Complex16> work(lwork);
	MKL_INT lrwork = MAX(1, 24 * n);
	Vec<double> rwork(lrwork);
	MKL_INT liwork = MAX(1, 10 * n);
	Vec<MKL_INT> iwork(liwork);
	zheevr(&jobz, &range, &uplo, &n, a, &n, &vl, &vu, &il, &iu, &abstol, &nr, w, z.p(), &n, isuppz.p(), 
		work.p(), &lwork, rwork.p(), &lrwork, iwork.p(), &liwork, &info);
	m = z.ct();
	return info;
}

#endif /* _MILAPACK_H_ */
