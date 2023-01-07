/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China) date 2013 - 2017
	Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2020 - 2022
*/

#ifndef _MIBLAS_H_
#define _MIBLAS_H_

#define MKL_Complex8  std::complex<float>
#define MKL_Complex16 std::complex<double>
#include <mkl.h>

// dot product lhs^+ dot rhs
inline double DOT(const Vec<double> &lhs, const Vec<double> &rhs)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(lhs.size(), rhs.size());
#endif
	double *X = lhs.p();
	double *Y = rhs.p();
	MKL_INT N = (MKL_INT)lhs.size();
	MKL_INT incX = 1;
	MKL_INT incY = 1;
	return cblas_ddot(N, X, incX, Y, incY);
}

// dot product lhs^+ dot rhs
inline MKL_Complex16 DOT(const Vec<MKL_Complex16> &lhs, const Vec<MKL_Complex16> &rhs)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(lhs.size(), rhs.size());
#endif
	MKL_Complex16 *X = lhs.p();
	MKL_Complex16 *Y = rhs.p();
	MKL_INT N = (MKL_INT)lhs.size();
	MKL_INT incX = 1;
	MKL_INT incY = 1;
	MKL_Complex16 dotc;
	cblas_zdotc_sub(N, (void *)X, incX, (void *)Y, incY, (void *)(&dotc));
	return dotc;
}

// matrix-vector multiplication y = a * x
inline void MUL(Vec<double> &y, const Mat<double> &a, const Vec<double> &x)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(y.size(), a.nrows(), a.ncols(), x.size());
#endif
	if (!y.size()) return;
	if (!x.size()) { y = 0.; return; }
	MKL_INT M = (MKL_INT)a.nrows();
	MKL_INT N = (MKL_INT)a.ncols();
	double *A = a.p();
	double *X = x.p();
	double *Y = y.p();
	double alpha = 1.;
	double beta = 0.;
	MKL_INT incX = 1;
	MKL_INT incY = 1;
	cblas_dgemv(CblasRowMajor, CblasNoTrans, M, N, alpha, A, N, X, incX, beta, Y, incY);
}
 
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// These two function was out of date cause OneAPI was come up, But It still have been used before 2022-4-19  !!!!!!!!!!!!!!
// (Waring Ifo:" warning #1478: function "mkl_dcsrmv";  warning #1478: function "mkl_dcsrgemv" "		      !!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//// sparse_matrix-vector multiplication y_{m,1} = a_{m,k} * x_{k,1}
//// y := alpha*A*x + beta*y
//// #not tested!
//inline void sparse_MUL(const Int& length_m, const Int& length_k, const Vec<Int>& roin, const Vec<Int>& colnu, const Vec<Real>& val, const Vec<Real>& KetVec, Vec<Real>& ret, const double alpha = 1., const double beta = 0.)
//{
//#ifdef _CHECK_DIMENSION_MATCH_
//	//ASSERT_EQ2(y.size(), a.nrows(), a.ncols(), x.size());
//	ASSERT_EQ(length_m + 1, roin.size());
//#endif
//
//	Vec<Int> pointerB(roin.truncate(0, length_m));
//	Vec<Int> piinterE(roin.truncate(1, length_m + 1));
//
//	MKL_INT* lm = (MKL_INT*)&length_m;
//	MKL_INT* lk = (MKL_INT*)&length_k;
//	Vec<MKL_INT> col(colnu);
//	Vec<MKL_INT> pntrb(pointerB);
//	Vec<MKL_INT> pntre(piinterE);
//
//	double* val_ = val.p();
//	MKL_INT* col_ = col.p();
//	MKL_INT* pntrb_ = pntrb.p();
//	MKL_INT* pntre_ = pntre.p();
//	double* x = KetVec.p();
//	double* y = ret.p();
//
//	mkl_dcsrmv("N", lm, lk, &alpha, "G", val_, col_, pntrb_, pntre_, x, &beta, y);
//}
//
//// sparse_matrix-vector multiplication y_{m,1} = a_{m,m} * x_{m,1}
//// y := A*x
//// tested!
//inline void sparse_MUL(Int dim, const Vec<Int>& roin, const Vec<Int>& colnu, const Vec<Real>& val, const Vec<Real>& KetVec, Vec<Real>& ret)
//{
//#ifdef _CHECK_DIMENSION_MATCH_
//	//ASSERT_EQ2(y.size(), a.nrows(), a.ncols(), x.size());
//	ASSERT_EQ(dim + 1, roin.size());
//#endif
//
//	MKL_INT* dimm = (MKL_INT*)&dim;
//	Vec<MKL_INT> col(colnu);
//	Vec<MKL_INT> rol(roin);
//
//	double* val_ = val.p();
//	MKL_INT* rol_ = rol.p();
//	MKL_INT* col_ = col.p();
//	double* x = KetVec.p();
//	double* y = ret.p();
//
//	mkl_dcsrgemv("N", dimm, val_, rol_, col_, x, y);
//}




// sparse_matrix-vector multiplication y_{m,1} = a_{m,k} * x_{k,1} 
// y := alpha*A*x + beta*y
inline void sparse_MUL(const Int& length_m, const Int& length_k, const Vec<Int>& roin, const Vec<Int>& colnu, const Vec<Real>& val, const Vec<Real>& KetVec, Vec<Real>& ret, const double alpha = 1., const double beta = 0.)
{
	Vec<Int> pointerB(roin.truncate(0, length_m));
	Vec<Int> piinterE(roin.truncate(1, length_m + 1));

	Vec<MKL_INT> col(colnu);
	Vec<MKL_INT> pntrb(pointerB);
	Vec<MKL_INT> pntre(piinterE);

	double* val_ = val.p();
	MKL_INT* col_ = col.p();
	MKL_INT* pntrb_ = pntrb.p();
	MKL_INT* pntre_ = pntre.p();

	struct matrix_descr descrA;
	sparse_matrix_t       csrA;

	double* x = KetVec.p();
	double* y = ret.p();

	mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ZERO, length_m, length_k, pntrb_, pntre_, col_, val_);
	descrA.type = SPARSE_MATRIX_TYPE_GENERAL;											// Create matrix descriptor
	// Analyze sparse matrix; choose proper kernels and workload balancing strategy
	mkl_sparse_optimize(csrA);
	// Compute y = alpha * A * x + beta * y
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, csrA, descrA, x, beta, y);	
	mkl_sparse_destroy(csrA);															// Release matrix handle and deallocate matrix
}


// C := op(A) *B
// The column indices of the output matrix (if in CSR format) can appear unsorted due to the algorithm chosen internally. To ensure sorted column indices (if that is important), call mkl_sparse_order().
inline void sparse_MUL(
	const Int& length_lr, const Int& length_lc, const Vec<Int>& roin_l, const Vec<Int>& colnu_l, const Vec<Real>& val_l,
	const Int& length_rr, const Int& length_rc, const Vec<Int>& roin_r, const Vec<Int>& colnu_r, const Vec<Real>& val_r,
													Vec<Int>& roin,			Vec<Int>& colnu,		 Vec<Real>& val)
{
	sparse_matrix_t	csrA = NULL, csrB = NULL, csrC = NULL;
	{// left side matrix data.

		double* val_ = val_l.p();
		MKL_INT* col_ = colnu_l.p();
		MKL_INT* pntrb_ = roin_l.p();
		MKL_INT* pntre_ = roin_l.p_one();
		mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ZERO, length_lr, length_lc, pntrb_, pntre_, col_, val_);
		mkl_sparse_optimize(csrA);				 // Analyze sparse matrix;
	}
	{// right side matrix data.

		double* val_ = val_r.p();
		MKL_INT* col_ = colnu_r.p();
		MKL_INT* pntrb_ = roin_r.p();
		MKL_INT* pntre_ = roin_r.p_one();
		mkl_sparse_d_create_csr(&csrB, SPARSE_INDEX_BASE_ZERO, length_rr, length_rc, pntrb_, pntre_, col_, val_);
		mkl_sparse_optimize(csrB);				 // Analyze sparse matrix;
	}

	mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, csrA, csrB, &csrC);
	// Release matrix handle and deallocate matrix
	mkl_sparse_destroy(csrA);
	mkl_sparse_destroy(csrB);
	// exort the sparse matrix C.
	sparse_index_base_t    indexing;
	MKL_INT  rows(length_lr), cols(length_rc);
	double* val_ = val.p();
	MKL_INT* col_ = colnu.p();
	MKL_INT* pntrb_ = roin.p();
	MKL_INT* pntre_ = roin.p_one();
	mkl_sparse_d_export_csr(csrC, &indexing, &rows, &cols, &pntrb_, &pntre_, &col_, &val_);

#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(rows, length_lr, cols, length_rc);
#endif
}

// C := alpha*op(A) + B
inline void sparse_ADD(
	const Int& length_lr, const Int& length_lc, const Vec<Int>& roin_l, const Vec<Int>& colnu_l, const Vec<Real>& val_l,
	const Int& length_rr, const Int& length_rc, const Vec<Int>& roin_r, const Vec<Int>& colnu_r, const Vec<Real>& val_r,
													  Vec<Int>& roin,		  Vec<Int>& colnu,		   Vec<Real>& val, 
	const double alpha = 1.)
{
#define CALL_AND_CHECK_STATUS(function, error_message) do { \
          if(function != SPARSE_STATUS_SUCCESS)             \
          {                                                 \
          printf(error_message); fflush(0);                 \
          status = 1;                                       \
          }                                                 \
} while(0)

	/* Declaration of values */
	double* values_A = NULL, * values_B = NULL, * values_C = NULL;
	MKL_INT* columns_A = NULL, * columns_B = NULL, * columns_C = NULL;
	MKL_INT* rowIndex_A = NULL, * rowIndex_B = NULL, * pointerB_C = NULL, * pointerE_C = NULL;

	sparse_matrix_t	csrA = NULL, csrB = NULL, csrC = NULL;
	MKL_INT status;
	// left side matrix data.

	Vec<MKL_INT> pntrb_a(roin_l.truncate(0, length_lr));
	Vec<MKL_INT> pntre_a(roin_l.truncate(1, length_lr + 1));
	Vec<MKL_INT> col_a(colnu_l);

	MKL_INT* pntrb_l = pntrb_a.p();
	MKL_INT* pntre_l = pntre_a.p();
	MKL_INT* col_l = col_a.p();
	double* va_l = val_l.p();
	mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ZERO, length_lr, length_lc, pntrb_l, pntre_l, col_l, va_l);
	CALL_AND_CHECK_STATUS(mkl_sparse_optimize(csrA), "Error after MKL_SPARSE_OPTIMIZE, csrA \n");
	// right side matrix data.

	Vec<MKL_INT> pntrb_b(roin_r.truncate(0, length_rr));
	Vec<MKL_INT> pntre_b(roin_r.truncate(1, length_rr + 1));
	Vec<MKL_INT> col_b(colnu_r);

	MKL_INT* pntrb_r = pntrb_b.p();
	MKL_INT* pntre_r = pntre_b.p();
	MKL_INT* col_r = col_b.p();
	double* va_r = val_r.p();
	mkl_sparse_d_create_csr(&csrB, SPARSE_INDEX_BASE_ZERO, length_rr, length_rc, pntrb_r, pntre_r, col_r, va_r);
	CALL_AND_CHECK_STATUS(mkl_sparse_optimize(csrB), "Error after MKL_SPARSE_OPTIMIZE, csrB \n");


	CALL_AND_CHECK_STATUS(mkl_sparse_d_add(SPARSE_OPERATION_NON_TRANSPOSE, csrA, alpha, csrB, &csrC), "Error after MKL_SPARSE_D_ADD  \n");

	CALL_AND_CHECK_STATUS(mkl_sparse_optimize(csrC), "Error after MKL_SPARSE_OPTIMIZE, csrC \n");
	// Release matrix handle and deallocate matrix
	mkl_sparse_destroy(csrA);
	mkl_sparse_destroy(csrB);
	// exort the sparse matrix C.
	sparse_index_base_t    indexing;
	MKL_INT  rows(length_lc), cols(length_rc);
	//DBG("Fine here.")
	CALL_AND_CHECK_STATUS(mkl_sparse_d_export_csr(csrC, &indexing, &rows, &cols, &pointerB_C, &pointerE_C, &columns_C, &values_C),	"Error after MKL_SPARSE_D_EXPORT_CSR  \n");
	Idx counter(0);
	for_Idx(i, 0, rows)for_Int(j, pointerB_C[i], pointerE_C[i])counter++;
	val.sm(counter, values_C);
	colnu.sm(counter, columns_C);
	roin.reset(rows + 1);
	for_Idx(i, 1, rows + 1)roin[i] = pointerE_C[i - 1];
	//DBG(NAV3(colnu,val,roin))
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(length_lr, length_rr, length_lc, length_rc);
#endif
	mkl_sparse_destroy(csrC);

}

// matrix-vector multiplication y = a * x
inline void MUL(Vec<MKL_Complex16> &y, const Mat<MKL_Complex16> &a, const Vec<MKL_Complex16> &x)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(y.size(), a.nrows(), a.ncols(), x.size());
#endif
	if (!y.size()) return;
	if (!x.size()) { y = 0.; return; }
	MKL_INT M = (MKL_INT)a.nrows();
	MKL_INT N = (MKL_INT)a.ncols();
	MKL_Complex16 *A = a.p();
	MKL_Complex16 *X = x.p();
	MKL_Complex16 *Y = y.p();
	MKL_Complex16 alpha = 1.;
	MKL_Complex16 beta = 0.;
	MKL_INT incX = 1;
	MKL_INT incY = 1;
	cblas_zgemv(CblasRowMajor, CblasNoTrans, 
		M, N, (void *)(&alpha), (void *)A, N, (void *)X, incX, (void *)(&beta), (void *)Y, incY);
}

// matrix-matrix multiplication a = b * c
inline void MUL(Mat<double> &a, const Mat<double> &b, const Mat<double> &c)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ3(a.nrows(), b.nrows(), b.ncols(), c.nrows(), c.ncols(), a.ncols());
#endif
	if (!a.nrows() || !a.ncols()) return;
	if (!b.ncols()) { a = 0.; return; }
	MKL_INT M = (MKL_INT)a.nrows();
	MKL_INT N = (MKL_INT)a.ncols();
	MKL_INT P = (MKL_INT)b.ncols();
	double alpha = 1.;
	double beta = 0.;
	double *A = b.p();
	double *B = c.p();
	double *C = a.p();
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
		M, N, P, alpha, A, P, B, N, beta, C, N);
}

// matrix-matrix multiplication a = b * c
inline void MUL(Mat<MKL_Complex16> &a, const Mat<MKL_Complex16> &b, const Mat<MKL_Complex16> &c)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ3(a.nrows(), b.nrows(), b.ncols(), c.nrows(), c.ncols(), a.ncols());
#endif
	if (!a.nrows() || !a.ncols()) return;
	if (!b.ncols()) { a = 0.; return; }
	MKL_INT M = (MKL_INT)a.nrows();
	MKL_INT N = (MKL_INT)a.ncols();
	MKL_INT P = (MKL_INT)b.ncols();
	MKL_Complex16 *A = b.p();
	MKL_Complex16 *B = c.p();
	MKL_Complex16 *C = a.p();
	MKL_Complex16 alpha = 1.;
	MKL_Complex16 beta = 0.;
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
		M, N, P, (void *)(&alpha), (void *)A, P, (void *)B, N, (void *)(&beta), (void *)C, N);
}

inline VecCmplx operator*(const MatReal &a, const VecCmplx &x) { return cmplx(a * real(x), a * imag(x)); }
inline VecCmplx operator*(const MatCmplx &a, const VecReal &x) { return a * cmplx(x); }


//	matrix-matrix multiplication a = alpha * b * c + beta * a
inline void MUL(MKL_Complex16 beta, Mat<MKL_Complex16> &a, MKL_Complex16 alpha, const Mat<MKL_Complex16> &b, const Mat<MKL_Complex16> &c)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ3(a.nrows(), b.nrows(), b.ncols(), c.nrows(), c.ncols(), a.ncols());
#endif
	if (!a.nrows() || !a.ncols()) return;
	if (!b.ncols()) { a = 0.; return; }
	MKL_INT M = (MKL_INT)a.nrows();
	MKL_INT N = (MKL_INT)a.ncols();
	MKL_INT P = (MKL_INT)b.ncols();
	MKL_Complex16 *A = b.p();
	MKL_Complex16 *B = c.p();
	MKL_Complex16 *C = a.p();
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		M, N, P, (void *)(&alpha), (void *)A, P, (void *)B, N, (void *)(&beta), (void *)C, N);
}

#endif /* _MIBLAS_H_ */
