/*
code developed by
	Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China) date 2013 - 2017
	Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2020 - 2022
*/

#ifndef _MAT_H_
#define _MAT_H_

#include "stdpfx.h"
#include "vec.h"

// class for matrix

// ========================= //
// see remarks for class Vec //
// ========================= //

template<typename T> inline void swapall(Mat<T> &lhs, Mat<T> &rhs);
template<typename T> inline std::ostream &operator<<(std::ostream &os, const Mat<T> &a);

template<typename T>
class Mat {
	friend void swapall<>(Mat<T> &lhs, Mat<T> &rhs);
    friend std::ostream &operator<< <>(std::ostream &os, const Mat<T> &a);
private:
	Idx m;				// # of rows
	Idx n;				// # of columns
	Int s;				// if borrow memory from elsewhere
	T* a;				// data, a[i * n + j] = Mat[i][j], consist with LAPACK_ROW_MAJOR
	Vec<Vec<T> > b;		// see Mat as Vec<Vec<T> >, only for use of subscript [i] and [i][j]
private:
	inline void refresh_b();								// refresh b when (m, n, a) changed
	inline void pdcam(Idx m_i, Idx n_i);					// pseudo delegating constructor, allocate memory
	inline void pdcbm(Idx m_i, Idx n_i, T *p);				// pseudo delegating constructor, borrow memory
public:
	typedef T value_type;									// make T available externally
	inline explicit Mat(): m(0), n(0), s(0), a(nullptr) {}	// default constructor
	inline explicit Mat(Idx n_i);							// set all elements with default, square matrix
	inline explicit Mat(Idx m_i, Idx n_i);					// set all elements with default
	inline explicit Mat(Idx m_i, Idx n_i, T t);				// set all elements with t
	// inline explicit Mat(Idx m_i, Idx n_i, T *p);			// borrow memory from p
	// inline explicit Mat(Idx m_i, Idx n_i, const Vec<T> &v);	// borrow memory from v, require m_i * n_i == v.size()
	inline Mat(const Mat &rhs);								// copy constructor
	Idx size() const { return m * n; }						// return matrix size, the number of elements
	Idx szof() const { return m * n * sizeof(T); }			// return matrix size in unit of char
	Idx nrows() const { return m; }
	Idx ncols() const { return n; }
	T* p() const { return a; }								// return initial address of matrix
	inline Vec<T> &operator[](Idx i);						// subscript to row i as a Vec borrowing memory from Mat
	inline const Vec<T> &operator[](Idx i) const;			// subscript to row i as a Vec borrowing memory from Mat
	inline Mat &operator=(const Mat &rhs);					// copy assignment operator, see the remarks above
	inline Mat &operator=(T t);								// assignment operator, Mat = t * I, require m == n
	inline Mat &reset();									// reset to status after invoke Mat()
	inline Mat &reset(Idx n_i);								// reset to status after invoke Mat(Idx n_i)
	inline Mat &reset(Idx m_i, Idx n_i);					// reset to status after invoke Mat(Idx m_i, Idx n_i)
	inline Mat &reset(Idx m_i, Idx n_i, T t);				// reset to status after invoke Mat(Idx m_i, Idx n_i, T t)
	inline Mat &reset(const Mat &rhs);						// reset to status after invoke Mat(const Mat &rhs)
	inline Mat &reset(const VEC<Vec<T>> &rhs);				// reset to status after invoke Mat(const Mat &rhs)
	inline Mat &sm(Idx m_i, Idx n_i, T *p);					// reset to status after invoke Mat(Idx m_i, Idx n_i, T *p), share memory with *p
	inline Mat &sm(Idx m_i, Idx n_i, Vec<T> &v);			// reset to status after invoke Mat(Idx m_i, Idx n_i, const Vec &v), share memory with v
	inline const Mat &sm(Idx m_i, Idx n_i, const T *p);		// reset to status after invoke Mat(Idx m_i, Idx n_i, T *p), share memory with *p
	inline const Mat &sm(Idx m_i, Idx n_i, const Vec<T> &v);// reset to status after invoke Mat(Idx m_i, Idx n_i, const Vec &v), share memory with v
	inline Mat &reshape(Idx m_i, Idx n_i);					// reshape, require m * n == m_i * n_i
	inline Mat &operator+=(const Mat &rhs);
	inline Mat &operator-=(const Mat &rhs);
	inline Mat &operator*=(T t);
	inline Mat &improve_unitarity();						// see note.LIOM.doc by Rong-Qiang He
	inline Mat operator-() const;
	inline Mat co() const;									// conjugate()
	inline Mat tr() const;									// transpose()
	inline Mat ct() const;									// conjugate_transpose()
	inline Mat truncate_row(Idx bgn, Idx end) const;
	inline Mat truncate(Idx bgn_i, Idx bgn_j, Idx end_i, Idx end_j) const;
	inline Vec<T> vec() const;								// convert a Vec
	inline Vec<T> diagonal() const;							// return diagonal elements as a Vec
	inline T trace() const;									// trace
	inline bool isunitary(Real eps = 100 * eps_Real) const;
	bool if_symmetric() const;
	inline Real norm_sqr() const { Vec<T> t; return t.sm(*this).norm_sqr(); }
	inline Real norm() const { Vec<T> t; return t.sm(*this).norm(); }
	inline Real norm_avg() const { return norm() / (SQRT(size()) + 1.E-128); }		// average norm per element
	inline Real max_abs_elem_diff(const Mat &b) const { return Vec<T>(*this).max_abs_elem_diff(Vec<T>(b)); }
	inline Real avg_abs_elem_diff(const Mat &b) const { return Vec<T>(*this).avg_abs_elem_diff(Vec<T>(b)); }
	inline ~Mat() { if (!s && a) delete [] a; }

	// deleting these two constructors so that the copy elision insecurity problem is avoided
	// inline explicit Mat(Idx m_i, Idx n_i, T *p);			// borrow memory from p
	// inline explicit Mat(Idx m_i, Idx n_i, const Vec<T> &v);	// borrow memory from v, require m_i * n_i == v.size()
};

// matrix types
typedef Mat<Int> MatInt;
typedef Mat<Idx> MatIdx;
typedef Mat<Real> MatReal;
typedef Mat<Cmplx> MatCmplx;

// related functions
template<typename T> inline Mat<T> dmat(Idx n, T d);
template<typename T> inline Mat<T> dmat(const Vec<T> &v);
template<typename T> inline void MUL(Vec<T> &y, const Mat<T> &a, const Vec<T> &x);
template<typename T> inline void MUL(Mat<T> &a, const Mat<T> &b, const Mat<T> &c);
template<typename T> inline T DOT(const Mat<T> &a, const Mat<T> &b);
template<typename T> inline T SUM(const Mat<T> &v);
template<typename T> inline T SUM_0toX(const Mat<T>& v,Idx row_bef, Int pstn);
template<typename T> inline T MIN(const Mat<T>& a) { Vec<T> t; return MIN(t.sm(a.size(), a.p())); };
template<typename T> inline T MAX(const Mat<T>& a) { Vec<T> t; return MAX(t.sm(a.size(), a.p())); };
template<typename T> inline Mat<T> SQRT(const Mat<T>& a);
template<typename T> inline Mat<T> operator*(const Mat<T> &lhs, const Mat<T> &rhs);
template<typename T> inline Vec<T> operator*(const Mat<T> &a, const Vec<T> &x);
template<typename T> inline Vec<T> operator*(const Vec<T> &x, const Mat<T> &a);
template<typename T> inline Mat<T> dmatmul(const Vec<T> &v, Mat<T> a);
template<typename T> inline Mat<T> dmatmul(Mat<T> a, const Vec<T> &v);
template<typename T> inline Mat<T> operator+(Mat<T> lhs, const Mat<T> &rhs) { return lhs += rhs; }
template<typename T> inline Mat<T> operator-(Mat<T> lhs, const Mat<T> &rhs) { return lhs -= rhs; }
template<typename T> inline Mat<T> operator*(Mat<T> z, T t) { return z *= t; }
template<typename T> inline Mat<T> operator*(T t, Mat<T> z) { return z *= t; }
template<typename T> inline bool operator==(const Mat<T> &lhs, const Mat<T> &rhs);
template<typename T> inline bool operator!=(const Mat<T> &lhs, const Mat<T> &rhs);
template<typename T> inline Mat<T> kronecker_product(const Mat<T> &a, const Mat<T> &b);
template<typename T> inline Mat<T> direct_sum(const Mat<T> &a, const Mat<T> &b);
template<typename T> inline Mat<T> commutator(const Mat<T> &a, const Mat<T> &b) { return a * b - b * a; }
template<typename T> inline Mat<T> anticommutator(const Mat<T> &a, const Mat<T> &b) { return a * b + b * a; }
template<typename T> inline Real RERR(const Mat<T> &a, const Mat<T> &b) { return (a - b).norm() * 2 / (a.norm() + b.norm() + 1.E-128); }
template<typename T> inline MatReal real(const Mat<T>& a);
template<typename T> inline MatReal imag(const Mat<T>& a);
inline MatCmplx cmplx(const MatReal& a, const MatReal& b);
inline MatCmplx cmplx(const MatReal& a);

// function definitions
template<typename T>
void swapall(Mat<T> &lhs, Mat<T> &rhs)
{
	SWAP(lhs.m, rhs.m);
	SWAP(lhs.n, rhs.n);
	SWAP(lhs.s, rhs.s);
	SWAP(lhs.a, rhs.a);
	swapall(lhs.b, rhs.b);
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const Mat<T> &a)
{
	os << std::endl;

	if (a.nrows() > 1024 || a.ncols() > 64) {
		return os << "matrix size = " << a.nrows() << " * " << a.ncols() << ", nrows > 1024 || ncols > 64" << std::endl;
	}

	const Int &wi = TOSTRLEN(a.nrows() - 1);
	const Int &wj = TOSTRLEN(a.ncols() - 1);
	Int w = 0;
	for_Idx (i, 0, a.nrows()) for_Idx (j, 0, a.ncols()) if (w < TOSTRLEN(a[i][j])) w = TOSTRLEN(a[i][j]);

	for_Idx (i, 0, a.nrows()) {
		os << "row[" << std::setw(wi) << i << "] =";
		for_Idx (j, 0, a.ncols()) os << " " << std::setw(w) << a[i][j];
		os << std::endl;
	}
	return os;
}

template<typename T>
// diagonal matrix
inline Mat<T> dmat(Idx n, T d)
{
	Mat<T> ret(n, n, T(0.));
	for_Idx (i, 0, n) ret[i][i] = d;
	return ret;
}

template<typename T>
// diagonal matrix
inline Mat<T> dmat(const Vec<T> &v)
{
	const Idx &n = v.size();
	Mat<T> ret(n, n, T(0.));
	for_Idx (i, 0, n) ret[i][i] = v[i];
	return ret;
}

template<typename T>
// pseudo delegating constructor, allocate memory
void Mat<T>::refresh_b()
{
	b.reset(m);
	for_Idx (i, 0, m) b[i].sm(n, a + i * n);
}

template<typename T>
// pseudo delegating constructor, allocate memory
void Mat<T>::pdcam(Idx m_i, Idx n_i)
{
	m = m_i;
	n = n_i;
	s = 0;
	a = new T[m * n];
	refresh_b();
}

template<typename T>
// pseudo delegating constructor, borrow memory
void Mat<T>::pdcbm(Idx m_i, Idx n_i, T *p)
{
	m = m_i;
	n = n_i;
	s = 1;
	a = p;
	refresh_b();
}

template<typename T>
// set all elements with default, square matrix
Mat<T>::Mat(Idx n_i)
{
	pdcam(n_i, n_i);
}

template<typename T>
// set all elements with default
Mat<T>::Mat(Idx m_i, Idx n_i)
{
	pdcam(m_i, n_i);
}

template<typename T>
// set all elements with t
Mat<T>::Mat(Idx m_i, Idx n_i, T t)
{
	pdcam(m_i, n_i);
	for_Idx (i, 0, m * n) a[i] = t;
}
/*
template<typename T>
// borrow memory from p
Mat<T>::Mat(Idx m_i, Idx n_i, T *p)
{
	pdcbm(m_i, n_i, p);
}

template<typename T>
// borrow memory from v, require m_i * n_i == v.size()
Mat<T>::Mat(Idx m_i, Idx n_i, const Vec<T> &v)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ1_PS(m_i * n_i, v.size(), NAVC2(m_i, n_i));
#endif
	pdcbm(m_i, n_i, v.p());
}
*/
template<typename T>
Mat<T>::Mat(const Mat &rhs)
{
	pdcam(rhs.m, rhs.n);
	for_Idx (i, 0, m * n) a[i] = rhs.a[i];
}

template<typename T>
Vec<T> &Mat<T>::operator[](Idx i)
{
#ifdef _CHECK_BOUNDS_
	if (i >= m) ERR("Mat row subscript out of bounds, subscript = " + STR(i) + ", m = " + STR(m) + ", n = " + STR(n));
#endif
	return b[i];
}

template<typename T>
const Vec<T> &Mat<T>::operator[](Idx i) const
{
#ifdef _CHECK_BOUNDS_
	if (i >= m) ERR("Mat row subscript out of bounds, subscript = " + STR(i) + ", m = " + STR(m) + ", n = " + STR(n));
#endif
	return b[i];
}

template<typename T>
Mat<T> &Mat<T>::operator=(const Mat &rhs)
{
	if (this == &rhs) return *this;
	if (!s && !a) pdcam(rhs.m, rhs.n);
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(m, rhs.m, n, rhs.n);
#endif
	for_Idx (i, 0, m * n) a[i] = rhs.a[i];
	return *this;
}

template<typename T>
// assignment operator, Mat = t * I, require m == n
Mat<T> &Mat<T>::operator=(const T t)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(m, n);
#endif
	for_Idx (i, 0, m * n) a[i] = T(0.);
	for_Idx (i, 0, m) a[i * n + i] = t;
	return *this;
}

template<typename T>
Mat<T> &Mat<T>::reset()
{
	if (!s && a) delete [] a;
	m = 0;
	n = 0;
	s = 0;
	a = nullptr;
	b.reset();
	return *this;
}

template<typename T>
Mat<T> &Mat<T>::reset(Idx n_i)
{
	if (!s && a) delete [] a;
	pdcam(n_i, n_i);
	return *this;
}

template<typename T>
Mat<T> &Mat<T>::reset(Idx m_i, Idx n_i)
{
	if (!s && a) delete [] a;
	pdcam(m_i, n_i);
	return *this;
}

template<typename T>
Mat<T> &Mat<T>::reset(Idx m_i, Idx n_i, T t)
{
	this->reset(m_i, n_i);
	for_Idx (i, 0, m * n) a[i] = t;
	return *this;
}

template<typename T>
Mat<T> &Mat<T>::sm(Idx m_i, Idx n_i, T *p)
{
	if (!s && a) delete [] a;
	pdcbm(m_i, n_i, p);
	return *this;
}

template<typename T>
const Mat<T> &Mat<T>::sm(Idx m_i, Idx n_i, const T *p)
{
	if (!s && a) delete [] a;
	pdcbm(m_i, n_i, p);
	return *this;
}

template<typename T>
Mat<T> &Mat<T>::sm(Idx m_i, Idx n_i, Vec<T> &v)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ1_PS(m_i * n_i, v.size(), NAVC2(m_i, n_i));
#endif
	this->sm(m_i, n_i, v.p());
	return *this;
}

template<typename T>
const Mat<T> &Mat<T>::sm(Idx m_i, Idx n_i, const Vec<T> &v)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ1_PS(m_i * n_i, v.size(), NAVC2(m_i, n_i));
#endif
	this->sm(m_i, n_i, v.p());
	return *this;
}

template<typename T>
Mat<T> &Mat<T>::reset(const Mat &rhs)
{
	this->reset();
	(*this) = rhs;
	return *this;
}

template<typename T>
Mat<T> &Mat<T>::reset(const VEC<Vec<T>> &rhs)
{
	this->reset(rhs.size(), rhs[0].size());
	for_Idx (i, 0, m)(*this)[i] = std::move(rhs[i]);
	return *this;
}

template<typename T>
Mat<T> &Mat<T>::reshape(Idx m_i, Idx n_i)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ1_PS(m * n, m_i * n_i, NAVC4(m, n, m_i, n_i));
#endif
	m = m_i;
	n = n_i;
	refresh_b();
	return *this;
}

template<typename T>
Mat<T> &Mat<T>::operator+=(const Mat &rhs)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(m, rhs.m, n, rhs.n);
#endif
	for_Idx (i, 0, m) for_Idx (j, 0, n) (*this)[i][j] += rhs[i][j];
	return *this;
}

template<typename T>
Mat<T> &Mat<T>::operator-=(const Mat &rhs)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(m, rhs.m, n, rhs.n);
#endif
	for_Idx (i, 0, m) for_Idx (j, 0, n) (*this)[i][j] -= rhs[i][j];
	return *this;
}

template<typename T>
// Mat *= t
Mat<T> &Mat<T>::operator*=(T t)
{
	for_Idx (i, 0, m * n) a[i] *= t;
	return *this;
}

template<typename T>
// Mat *= t
Mat<T> &Mat<T>::improve_unitarity()
{
	Mat<T> &U = *this;
	return U = 1.5 * U - 0.5 * U * U.ct() * U;
}

template<typename T>
Mat<T> Mat<T>::operator-() const
{
	Mat<T> ret(m, n);
	for_Idx (i, 0, m) for_Idx (j, 0, n) ret[i][j] = -(*this)[i][j];
	return ret;
}

template<typename T>
Mat<T> Mat<T>::co() const
{
	Mat<T> ret(m, n);
	for_Idx (i, 0, m) for_Idx (j, 0, n) ret[i][j] = cnjg((*this)[i][j]);
	return ret;
}

template<typename T>
Mat<T> Mat<T>::tr() const
{
	Mat<T> ret(n, m);
	for_Idx (i, 0, m) for_Idx (j, 0, n) ret[j][i] = (*this)[i][j];
	return ret;
}

template<typename T>
Mat<T> Mat<T>::ct() const
{
	Mat<T> ret(n, m);
	for_Idx (i, 0, m) for_Idx (j, 0, n) ret[j][i] = cnjg((*this)[i][j]);
	return ret;
}

template<typename T>
Mat<T> Mat<T>::truncate_row(Idx bgn, Idx end) const
{
	if (bgn < end) {
		Mat<T> b(end - bgn, this->ncols());
		for_Idx(i, 0, b.nrows()) {
			for_Idx(j, 0, b.ncols()) {
				b[i][j] = (*this)[bgn + i][j];
			}
		}
		return b;
	}
	else {
		return Mat<T>();
	}
}

template<typename T>
Mat<T> Mat<T>::truncate(Idx bgn_i, Idx bgn_j, Idx end_i, Idx end_j) const
{
	if (bgn_i < end_i && bgn_j < end_j) {
		Mat<T> b(end_i - bgn_i, end_j - bgn_j);
		for_Idx(i, 0, b.nrows()) {
			for_Idx(j, 0, b.ncols()) {
				b[i][j] = (*this)[bgn_i + i][bgn_j + j];
			}
		}
		return b;
	}
	else {
		return Mat<T>();
	}
}

template<typename T>
// convert to a Vec
Vec<T> Mat<T>::vec() const
{
	Vec<T> v(m * n);
	Vec<T> t;
	v = t.sm(*this);
	return v;
}

template<typename T>
Vec<T> Mat<T>::diagonal() const
{
	Vec<T> v(MIN(m, n));
	for_Idx (i, 0, v.size()) v[i] = (*this)[i][i];
	return v;
}

template<typename T>
T Mat<T>::trace() const
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(m, n);
#endif
	T sum = T(0.);
	for_Idx (i, 0, m) sum += (*this)[i][i];
	return sum;
}

template<typename T>
bool Mat<T>::isunitary(Real eps) const
{
	Mat<T> a = m > n ? (*this).ct() * (*this) : (*this) * (*this).ct();
	Real thh = SQRT(a.nrows()) * eps;
	for_Idx (i, 0, a.nrows()) a[i][i] -= T(1.);
	Vec<T> v(a);
	for_Idx (i, 0, v.size()) if (ABS(v[i]) > thh) return false;
	return true;
}
template<typename T>
inline bool Mat<T>::if_symmetric() const
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ((*this).ncols(), (*this).nrows()); // on support for one phalanx matrix.
#endif
	for_Idx(i, 0, m) for_Idx(j, 0, n) if (ABS(((*this)[j][i] - (*this)[i][j])) > 1.E-9) {
		WRN("This position were not equal:"+NAV2((*this)[j][i], (*this)[i][j]))
		return false;
	}
	return true;
}
/*
// use MKL instead
template<typename T>
//	matrix-vector multiplication y = a * x
inline void MUL(Vec<T> &y, const Mat<T> &a, const Vec<T> &x)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(y.size(), a.nrows(), a.ncols(), x.size());
#endif
	const Idx &m = a.nrows();
	const Idx &n = a.ncols();
	for_Idx (i, 0, m) { T t = T(0.); for_Idx (j, 0, n) t += a[i][j] * x[j]; y[i] = t; }
}
// use MKL instead
template<typename T>
//	matrix-matrix multiplication a = b * c
inline void MUL(Mat<T> &a, const Mat<T> &b, const Mat<T> &c)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ3(a.nrows(), b.nrows(), b.ncols(), c.nrows(), c.ncols(), a.ncols());
#endif
	const Idx &l = a.nrows();
	const Idx &m = b.ncols();
	const Idx &n = a.ncols();
	for_Idx (i, 0, l) for_Idx (k, 0, n) { T t = T(0.); for_Idx (j, 0, m) t += b[i][j] * c[j][k]; a[i][k] = t; }
}
*/
template<typename T>
inline T DOT(const Mat<T> &a, const Mat<T> &b)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(a.nrows(), b.nrows(), a.ncols(), b.ncols());
#endif
	return DOT(Vec<T>(a), Vec<T>(b));
}


template<typename T>
inline Mat<T> SQRT(const Mat<T>& a)
{
	Mat<T> b(a);
	for_Idx(i, 0, b.nrows()) for_Idx(j, 0, b.ncols()) {
		b[i][j] = SQRT(a[i][j]);
	}
	return b;
}

template<typename T>
inline Mat<T> operator*(const Mat<T> &lhs, const Mat<T> &rhs)
{
	Mat<T> ret(lhs.nrows(), rhs.ncols());
	MUL(ret, lhs, rhs);
	return ret;
}

template<typename T>
inline Vec<T> operator*(const Mat<T> &a, const Vec<T> &x)
{
	Vec<T> ret(a.nrows());
	MUL(ret, a, x);
	return ret;
}

template<typename T>
inline Vec<T> operator*(const Vec<T> &x, const Mat<T> &a)
{
	Vec<T> ret(a.nrows());
	MUL(ret, a.tr(), x);
	return ret;
}

template<typename T>
inline Mat<T> dmatmul(const Vec<T> &v, Mat<T> a)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(v.size(), a.nrows());
#endif
	for_Idx (i, 0, a.nrows()) {
		a[i] *= v[i];
	}
	return a;
}

template<typename T>
inline Mat<T> dmatmul(Mat<T> a, const Vec<T> &v)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(v.size(), a.ncols());
#endif
	for_Idx (i, 0, a.nrows()) {
		a[i] *= v;
	}
	return a;
}

template<typename T>
inline bool operator==(const Mat<T> &lhs, const Mat<T> &rhs)
{
	if (lhs.nrows() != rhs.nrows()) return 0;
	if (lhs.ncols() != rhs.ncols()) return 0;
	for_Idx (i, 0, lhs.nrows()) for_Idx (j, 0, lhs.ncols()) if (lhs[i][j] != rhs[i][j]) return 0;
	return 1;
}

template<typename T>
inline bool operator!=(const Mat<T> &lhs, const Mat<T> &rhs)
{
	return !(lhs == rhs);
}

template<typename T>
inline Mat<T> kronecker_product(const Mat<T> &a, const Mat<T> &b)
{
	Mat<T> ret(a.nrows() * b.nrows(), a.ncols() * b.ncols());
	for_Idx (i, 0, a.nrows()) for_Idx (k, 0, b.nrows()) {
		Idx ik = i * b.nrows() + k;
		for_Idx (j, 0, a.ncols()) for_Idx (l, 0, b.ncols()) {
			Idx jl = j * b.ncols() + l;
			ret[ik][jl] = a[i][j] * b[k][l];
		}
	}
	return ret;
}

template<typename T>
inline Mat<T> direct_sum(const Mat<T> &a, const Mat<T> &b)
{
    Mat<T> ret(a.nrows() + b.nrows(), a.ncols() + b.ncols(), 0.);
    for_Idx (i, 0, a.nrows()) for_Idx (j, 0, a.ncols()) {
        ret[i][j] = a[i][j];
    }    
	for_Idx (i, 0, b.nrows()) for_Idx (j, 0, b.ncols()) {
        ret[i + a.nrows()][j + a.ncols()] = b[i][j];
    }
	return ret;
}


template<typename T>
inline T SUM(const Mat<T> &v)
{
	T ret = T(0.);
	for_Idx(i, 0, v.nrows()) {
		for_Idx(j, 0, v.ncols()) {
			ret += v[i][j];
		}
	}
	return ret;
}

template<typename T>
inline T SUM_0toX(const Mat<T>& v,Idx row_bef, Int pstn)
{
	ASSERT_GE(v.size(), pstn);
	T ret = T(0.);
	for_Idx(i, 0, row_bef) ret += SUM(v[i]);
	for_Idx(j, 0, pstn) ret += v[row_bef][j];
	return ret;
}

template<typename T>
inline MatReal real(const Mat<T>& a)
{
	MatReal r(a.nrows(), a.ncols());
	for_Idx(i, 0, r.nrows()) {
		for_Idx(j, 0, r.ncols()) {
			r[i][j] = real(a[i][j]);
		}
	}
	return r;
}

template<typename T>
inline MatReal imag(const Mat<T>& a)
{
	MatReal r(a.nrows(), a.ncols());
	for_Idx(i, 0, r.nrows()) {
		for_Idx(j, 0, r.ncols()) {
			r[i][j] = imag(a[i][j]);
		}
	}
	return r;
}

inline MatCmplx cmplx(const MatReal& a, const MatReal& b)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(a.nrows(), b.nrows(), a.ncols(), b.ncols());
#endif
	MatCmplx r(a.nrows(), a.ncols());
	for_Idx(i, 0, r.nrows()) {
		for_Idx(j, 0, r.ncols()) {
			r[i][j] = cmplx(a[i][j], b[i][j]);
		}
	}
	return r;
}

inline MatCmplx cmplx(const MatReal& a)
{
	MatCmplx r(a.nrows(), a.ncols());
	for_Idx(i, 0, r.nrows()) {
		for_Idx(j, 0, r.ncols()) {
			r[i][j] = cmplx(a[i][j], 0.);
		}
	}
	return r;
}

#endif /* _MAT_H_ */
