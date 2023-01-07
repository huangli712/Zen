/*
code developed by
	Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China) date 2013 - 2017
	Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2020 - 2022
*/

#ifndef _VEC_H_
#define _VEC_H_

#include "stdpfx.h"

// class for vector

// remarks for class Vec, which apply similarly to class Mat

// ============================================ Caveat ============================================ //
// Copy elision is the only allowed form of optimization that can change the observable side-effects. 
// Because some compilers do not perform copy elision in every situation where it is allowed (e.g., in debug mode), 
// programs that rely on the side-effects of copy/move constructors and destructors are not portable.
// Even when copy elision takes place and the copy-/move-constructor is not called, it must be present 
// and accessible (as if no optimization happened at all), otherwise the program is ill-formed.
// http://en.cppreference.com/w/cpp/language/copy_elision
// ================================================================================================ //

// this means that, for example, VecReal v = VecReal(mat) will be changed into
// VecReal v(mat) after compiler optimization, where mat is an object of MatReal. 
// This dangerous operation can be avoided by using VecReal v; v = VecReal(mat);

// ================================================================================================ //
// the constructors and the copy assignment operator are very unusual, please study them carefully! //
// ================================================================================================ //


// a Vec may borrow memory from outside, another Vec, a Mat, an array, or a pointer, etc.
// this is necessary because we often want to treat a consective part of a big object as a Vec.
// we should not establish another vector class that has its own resource.


// copy assignment operator is call when assigning an object to an existing object
// an on-state vector can not be assigned with a vector of different size
// off-state <=> [state after Vec() or reset()] <=> (s == 0 && a == nullptr) and <NOT=> any other state
// so we can use reset() to turn an on-state vector into an off-state vector


// copy constructor is called when passed or returned by value, or constructed by another object of the same class.
// in brief, copy constructor is called when a new object is created with an existing object.
// please preferentially make Vec parameters be passed or returned by reference rather than value


// copy assignment operator might be replaced by copy constructor when compiler optimizes
// this is why we place explicit specifier in front of constructors. but this problem is not clear!
// in C++11 multi-parameter constructors can be implicitly converted to with brace initialization.
// however, before C++11 explicit only applied to single-argument constructors.
// for multiple-argument constructors, it was ignored and had no effect.
// http://stackoverflow.com/questions/1118680/explicit-keyword-on-multi-arg-constructor
// http://stackoverflow.com/questions/4467142/why-is-explicit-allowed-for-default-constructors-and-constructors-with-2-or-more/4467658#4467658
// https://msdn.microsoft.com/en-us/library/wwywka61.aspx

// C++ note
// I find that implcit template class instantiations and member function instantiations
// can call functions after template class and member function definitions. An example is that
// functions DOT and MUL using MKL can be called by Vec and Mat memeber functions.
// This is dangerous, so we want to avoid this.


template<typename T> class Vec;
template<typename T> class Mat;
template<typename T> inline void swapall(Vec<T> &lhs, Vec<T> &rhs);
template<typename T> inline std::ostream &operator<<(std::ostream &os, const Vec<T> &v);

template<typename T>
class Vec {
	friend void swapall<>(Vec<T> &lhs, Vec<T> &rhs);
    friend std::ostream &operator<< <>(std::ostream &os, const Vec<T> &v);
private:
	Idx n;		// size of array
	Int s;		// if borrow memory from elsewhere
	T *a;		// data, a[i], i is in [0, n)
public:
	typedef T value_type;								// make T available externally
	typedef VEC<T> StdVec;								// for convenience
	explicit Vec(): n(0), s(0), a(nullptr) {}
	explicit Vec(Idx n_i): n(n_i), s(0), a(new T[n]) {}													// set all elements with default
	explicit Vec(Idx n_i, T t): n(n_i), s(0), a(new T[n]) { for_Idx (i, 0, n) a[i] = t; }				// set all elements with t
	// explicit Vec(Idx n_i, T *p): n(n_i), s(1), a(p) {}												// borrow memory from p
	// explicit Vec(const Mat<T> &m): n(m.size()), s(1), a(m.p()) {}									// borrow memory from Mat
	explicit Vec(const StdVec &v): n(v.size()), s(0), a(new T[n]) { for_Idx (i, 0, n) a[i] = v[i]; }	// convert VEC into Vec
	explicit Vec(const std::initializer_list<T>& v) : n(v.size()), s(0), a(new T[n]) {
		Idx count = 0;
		for (auto item : v) a[count++] = item;
	}
	Vec(const Vec &v): n(v.n), s(0), a(new T[n]) { for_Idx (i, 0, n) a[i] = v[i]; }						// copy constructor
	StdVec stdvec() const { return StdVec(a, a + n); }	// convert Vec into VEC
	Idx size() const { return n; }						// return vector size
	Idx szof() const { return n * sizeof(T); }			// return vector size in unit of char
	T *p() const { return a; }							// return initial address of vector
	T *p_one() const { return a + 1; }					// return (initial + 1) address of vector
	T *begin() const { return a; }						// v.begin()
	T *end() const { return a + n; }					// v.end()
	inline T &operator[](Idx i);						// subscript operator
	inline const T &operator[](Idx i) const;			// subscript operator
	inline Vec &operator=(const Vec &rhs);				// copy assignment operator, see the remarks above
	inline Vec &operator=(const StdVec &rhs);			// assign VEC to Vec
	inline Vec &operator=(T t);							// assignment operator, a[i] = t for all i
	inline Vec &reset();								// reset to status after invoke Vec()
	inline Vec &reset(Idx n_i);							// reset to status after invoke Vec(Idx n_i)
	inline Vec &reset(Idx n_i, T t);					// reset to status after invoke Vec(Idx n_i, T t)
	inline Vec &reset(const StdVec &v);					// reset to status after invoke Vec(const StdVec &v)
	inline Vec &reset(const Vec &v);					// reset to status after invoke Vec(const Vec &rhs)
	inline Vec &sm(Idx n_i, T *p);						// reset to status after invoke Vec(Idx n_i, T *p), share memory with *p
	inline Vec &sm(Mat<T> &x);							// reset to status after invoke Vec(const Mat<T> &x), share memory with x
	inline const Vec &sm(Idx n_i, const T *p);			// reset to status after invoke Vec(Idx n_i, T *p), share memory with *p
	inline const Vec &sm(const Mat<T> &x);				// reset to status after invoke Vec(const Mat<T> &x), share memory with x
	inline Vec &operator+=(const Vec &rhs);
	inline Vec &operator-=(const Vec &rhs);
	inline Vec &operator*=(const Vec &rhs);
	inline Vec &operator*=(const Mat<T> &rhs);
	inline Vec &operator+=(T t);
	inline Vec &operator-=(T t);
	inline Vec &operator*=(T t);
	inline Vec &operator*(const Vec<Bool>& rhs);
	inline Vec operator-() const;
	inline Idx idx_min() const;
	inline Idx idx_max() const;
	inline std::string string() const;
	inline Vec truncate(Idx bgn, Idx end) const;
	inline Vec reverse() const;
	inline Vec co() const;								// conjugate()
	inline Mat<T> mat(Idx m_i, Idx n_i) const;			// convert to a Mat
	inline Real norm_sqr() const { return real(DOT(*this, *this)); }
	inline Real norm() const { return SQRT(norm_sqr()); }
	inline Real norm_avg() const { return norm() / (SQRT(size()) + 1.E-128); }		// average norm per element
	inline Vec& normalize() { return (*this) *= T(INV(norm())); }
	inline bool isnormalized() const { return ABS(norm() - 1.) < 100 * SQRT(n) * eps_Real; }
	inline Real max_abs_elem_diff(const Vec &b) const { return MAX(ABS(b - *this)); }
	inline Real avg_abs_elem_diff(const Vec &b) const { return AVG(ABS(b - *this)); }
	inline ~Vec() { if (!s && a) delete [] a; }

	// deleting these two constructors so that the copy elision insecurity problem is avoided
	// explicit Vec(Idx n_i, T *p): n(n_i), s(1), a(p) {}													// borrow memory from p
	// explicit Vec(const Mat<T> &m): n(m.size()), s(1), a(m.p()) {}										// borrow memory from Mat
};

// vector types
typedef Vec<Int> VecInt;
typedef Vec<Idx> VecIdx;
typedef Vec<Real> VecReal;
typedef Vec<Cmplx> VecCmplx;
typedef Vec<bool> VecBool;
typedef std::vector<Int> StdVecInt;
typedef std::vector<Idx> StdVecIdx;
typedef std::vector<Real> StdVecReal;
typedef std::vector<Cmplx> StdVecCmplx;
typedef std::vector<Str> StdVecStr;

// related functions
template<typename T> inline Vec<T> operator+(Vec<T> lhs, const Vec<T> &rhs) { return lhs += rhs; }
template<typename T> inline Vec<T> operator-(Vec<T> lhs, const Vec<T> &rhs) { return lhs -= rhs; }
template<typename T> inline Vec<T> operator*(Vec<T> lhs, const Vec<T> &rhs) { return lhs *= rhs; }
template<typename T> inline Vec<T> operator*(Vec<T> lhs, T t) { return lhs *= t; }
template<typename T> inline Vec<T> operator*(T t, Vec<T> rhs) { return rhs *= t; }
template<typename T> inline Vec<T> concat(const Vec<T> &lhs, const Vec<T> &rhs);
template<typename T> inline bool operator==(const Vec<T> &lhs, const Vec<T> &rhs);
template<typename T> inline bool operator!=(const Vec<T> &lhs, const Vec<T> &rhs);
template<typename T> inline Mat<T> outer_product(const Vec<T> &lhs, const Vec<T> &rhs);
template<typename T> inline T DOT(const Vec<T> &lhs, const Vec<T> &rhs);
template<typename T> inline T SUM(const Vec<T> &v);
template<typename T> inline T SUM_0toX(const Vec<T>& v, Int pstn);
template<typename T> inline Vec<T> VECVectoVec(const VEC<Vec<T>>& m);
template<typename T> inline T MIN(const Vec<T> &v);
template<typename T> inline T MAX(const Vec<T> &v);
template<typename T> inline T AVG(const Vec<T> &v);
template<typename T> inline T CHI(const Vec<T> &v);
template<typename T> inline T DEV(const Vec<T> &v);
template<typename T> inline Real MAG(const Vec<T> &v);
template<typename T> inline Vec<T> ABS(const Vec<T> &a);
inline VecReal ABS(const VecCmplx &a);
template<typename T> inline Vec<T> EXP(const Vec<T> &a);
template<typename T> inline Vec<T> SQR(const Vec<T> &a);
template<typename T> inline Vec<T> SQRT(const Vec<T> &a);
template<typename T> inline Vec<T> INV(const Vec<T> &a);
template<typename T> inline Vec<T> kronecker_product(const Vec<T> &a, const Vec<T> &b);
template<typename T> inline Real RERR(const Vec<T> &a, const Vec<T> &b) { return (a - b).norm() * 2 / (a.norm() + b.norm() + 1.E-128); }
template<typename T> inline VecReal real(const Vec<T> &v);
template<typename T> inline VecReal imag(const Vec<T> &v);
inline VecCmplx cmplx(const VecReal &a, const VecReal &b);
inline VecCmplx cmplx(const VecReal &a);

// function definitions
template<typename T>
void swapall(Vec<T> &lhs, Vec<T> &rhs)
{
	swapall(lhs.n, rhs.n);
	swapall(lhs.s, rhs.s);
	swapall(lhs.a, rhs.a);
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const Vec<T> &v)
{
	os << std::endl;
	if (v.size() > 1024) return os << "vector size = " << v.size() << " > 1024" << std::endl;
	const Int &wi = TOSTRLEN(v.size() - 1);
	Int w = 0;
	for_Idx (i, 0, v.size()) if (w < TOSTRLEN(v[i])) w = TOSTRLEN(v[i]);
	for_Idx (i, 0, v.size()) os << "[" << std::setw(wi) << i << "] = " << std::setw(w) << v[i] << std::endl;
	return os;
}

template<typename T>
inline T &Vec<T>::operator[](Idx i)
{
#ifdef _CHECK_BOUNDS_
	if (i >= n) ERR("Vec subscript out of bounds, subscript = " + STR(i) + ", vector size = " + STR(n));
#endif
	return a[i];
}

template<typename T>
inline const T &Vec<T>::operator[](Idx i) const
{
#ifdef _CHECK_BOUNDS_
	if (i >= n) ERR("Vec subscript out of bounds, subscript = " + STR(i) + ", vector size = " + STR(n));
#endif
	return a[i];
}

template<typename T>
Vec<T> &Vec<T>::operator=(const Vec &rhs)
{
	if (this == &rhs) return *this;
	if (!rhs.s && !rhs.a) return (*this).reset();
	if (!s && !a) { n = rhs.n; a = new T[n]; }
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(n, rhs.n);
#endif
	for_Idx (i, 0, n) a[i] = rhs[i];
	return *this;
}

template<typename T>
Vec<T> &Vec<T>::operator=(const StdVec &rhs)
{
	if (!s && !a) { n = rhs.size(); a = new T[n]; }
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(n, rhs.size());
#endif
	for_Idx (i, 0, n) a[i] = rhs[i];
	return *this;
}

template<typename T>
Vec<T> &Vec<T>::operator=(T t)
{
	for_Idx (i, 0, n) a[i] = t;
	return *this;
}

template<class T>
Vec<T> &Vec<T>::reset()
{
	if (!s && a) delete [] a;
	n = 0;
	s = 0;
	a = nullptr;
	return *this;
}

template<class T>
Vec<T> &Vec<T>::reset(Idx n_i)
{
	if (!s && a) delete [] a;
	n = n_i;
	s = 0;
	a = new T[n];
	return *this;
}

template<typename T>
Vec<T> &Vec<T>::reset(Idx n_i, T t)
{
	this->reset(n_i);
	for_Idx (i, 0, n) a[i] = t;
	return *this;
}

template<typename T>
Vec<T> &Vec<T>::sm(Idx n_i, T *p)
{
	if (!s && a) delete [] a;
	n = n_i;
	s = 1;
	a = p;
	return *this;
}

template<typename T>
const Vec<T> &Vec<T>::sm(Idx n_i, const T *p)
{
	if (!s && a) delete [] a;
	n = n_i;
	s = 1;
	a = p;
	return *this;
}

template<typename T>
Vec<T> &Vec<T>::sm(Mat<T> &x)
{
	this->sm(x.size(), x.p());
	return *this;
}

template<typename T>
const Vec<T> &Vec<T>::sm(const Mat<T> &x)
{
	this->sm(x.size(), x.p());
	return  *this;
}

template<typename T>
Vec<T> &Vec<T>::reset(const StdVec &v)
{
	this->reset();
	*this = v;
	return *this;
}

template<typename T>
Vec<T> &Vec<T>::reset(const Vec &v)
{
	this->reset();
	*this = v;
	return *this;
}

template<typename T>
Vec<T> &Vec<T>::operator+=(const Vec &rhs)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(n, rhs.n);
#endif
	for_Idx (i, 0, n) (*this)[i] += rhs[i];
	return *this;
}

template<typename T>
Vec<T> &Vec<T>::operator-=(const Vec &rhs)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(n, rhs.n);
#endif
	for_Idx (i, 0, n) (*this)[i] -= rhs[i];
	return *this;
}

template<typename T>
Vec<T> &Vec<T>::operator*=(const Vec &rhs)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(n, rhs.n);
#endif
	for_Idx (i, 0, n) (*this)[i] *= rhs[i];
	return *this;
}

template<typename T>
Vec<T> &Vec<T>::operator*=(const Mat<T> &rhs)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(n, rhs.nrows());
#endif
	Vec<T> tmp = *this;
	MUL(*this, rhs.tr(), tmp);
	return *this;
}

template<typename T>
Vec<T>& Vec<T>::operator*(const Vec<Bool>& rhs)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(n, rhs.size());
#endif
	for_Idx(i, 0, n)(*this)[i] = (rhs[i] ? (*this)[i] : 0.);
	return (*this);
}

template<typename T> Vec<T> &Vec<T>::operator+=(T t) { for_Idx (i, 0, n) (*this)[i] += t; return *this; }
template<typename T> Vec<T> &Vec<T>::operator-=(T t) { for_Idx (i, 0, n) (*this)[i] -= t; return *this; }
template<typename T> Vec<T> &Vec<T>::operator*=(T t) { for_Idx (i, 0, n) (*this)[i] *= t; return *this; }

template<typename T>
Vec<T> Vec<T>::operator-() const
{
	Vec<T> ret(n);
	for_Idx (i, 0, n) ret[i] = -a[i];
	return ret;
}

template<typename T>
Idx Vec<T>::idx_min() const
{
	ASSERT_NE(n, 0);
	T tmp = a[0];
	Idx idx = 0;
	for_Idx(i, 1, n) if (tmp > a[i]) {
		tmp = a[i];
		idx = i;
	}
	return idx;
}

template<typename T>
Idx Vec<T>::idx_max() const
{
	ASSERT_NE(n, 0);
	T tmp = a[0];
	Idx idx = 0;
	for_Idx(i, 1, n) if (tmp < a[i]) {
		tmp = a[i];
		idx = i;
	}
	return idx;
}

template<typename T>
std::string Vec<T>::string() const
{
	std::string ret;
	for_Idx(i, 0, n) ret += char(a[i]);
	return ret;
}

template<typename T>
Vec<T> Vec<T>::truncate(Idx bgn, Idx end) const
{
	if (bgn < end) {
		Vec<T> b(end - bgn);
		for_Idx(i, 0, b.size()) {
			b[i] = (*this)[bgn + i];
		}
		return b;
	}
	else {
		return Vec<T>();
	}
}

template<typename T>
Vec<T> Vec<T>::reverse() const
{
	Vec<T> ret(n);
	for_Idx (i, 0, n) ret[i] = a[n - i - 1];
	return ret;
}

template<typename T>
Vec<T> Vec<T>::co() const
{
	Vec<T> ret(n);
	for_Idx (i, 0, n) ret[i] = cnjg(a[i]);
	return ret;
}

template<typename T>
Mat<T> Vec<T>::mat(Idx m_i, Idx n_i) const
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ1_PS(m_i * n_i, n, NAVC2(m_i, n_i));
#endif
	Mat<T> a(m_i, n_i);
	Vec<T> v;
	v.sm(a) = *this;
	return a;
}

template<typename T>
inline Vec<T> concat(const Vec<T> &lhs, const Vec<T> &rhs)
{
	Vec<T> ret(lhs.size() + rhs.size());
	Vec<T> t;
	t.sm(lhs.size(), ret.p()) = lhs;
	t.sm(rhs.size(), ret.p() + lhs.size()) = rhs;
	return ret;
}

template<typename T>
inline bool operator==(const Vec<T> &lhs, const Vec<T> &rhs)
{
	if (lhs.size() != rhs.size()) return 0;
	for_Idx (i, 0, lhs.size()) if (lhs[i] != rhs[i]) return 0;
	return 1;
}

template<typename T>
inline bool operator!=(const Vec<T> &lhs, const Vec<T> &rhs)
{
	return !(lhs == rhs);
}

template<typename T>
// lhs * rhs^dagger
inline Mat<T> outer_product(const Vec<T> &lhs, const Vec<T> &rhs)
{
	Mat<T> ret(lhs.size(), rhs.size());
	for_Idx (i, 0, lhs.size()) for_Idx (j, 0, rhs.size()) ret[i][j] = lhs[i] * cnjg(rhs[j]);
	return ret;
}

/*
template<typename T>
inline T DOT(const Vec<T> &lhs, const Vec<T> &rhs)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(lhs.size(), rhs.size());
#endif
	T ret = 0;
	for_Idx (i, 0, lhs.size()) ret += cnjg(lhs[i]) * rhs[i];
	return ret;
}
*/

template<typename T>
inline T SUM(const Vec<T> &v)
{
	T ret = T(0.);
	for_Idx (i, 0, v.size()) ret += v[i];
	return ret;
}

template<typename T>
inline T SUM_0toX(const Vec<T>& v, Int pstn)
{
	ASSERT_GE(v.size(), pstn);
	T ret = T(0.);
	for_Idx(i, 0, pstn) ret += v[i];
	return ret;
}

template<typename T>
inline Vec<T> VECVectoVec(const VEC<Vec<T>>& m)
{
	VEC<T> a;
	for (const auto& mvec_i : m) {
		for_Int(i, 0, mvec_i.size()) {
			a.push_back(mvec_i[i]);
		}
	}
	Vec<T> ret(a);
	return ret;
}

template<typename T>
inline T MIN(const Vec<T> &v)
{
	ASSERT_NE(v.size(), 0);
	T ret = v[0];
	for_Idx (i, 1, v.size()) if (ret > v[i]) ret = v[i];
	return ret;
}

template<typename T>
inline T MAX(const Vec<T> &v)
{
	ASSERT_NE(v.size(), 0);
	T ret = v[0];
	for_Idx (i, 1, v.size()) if (ret < v[i]) ret = v[i];
	return ret;
}

template<typename T>
inline T AVG(const Vec<T> &v)
{
	return v.size() > 0 ? SUM(v) / v.size() : T(0.);
}

template<typename T>
inline T CHI(const Vec<T> &v)
{
	T avg = AVG(v);
	T sum = T(0.);
	for_Idx (i, 0, v.size()) sum += cnjg(v[i] - avg) * (v[i] - avg);
	return sum;
}

template<typename T> inline T DEV(const Vec<T> &v)
{
	return v.size() > 0 ? sqrt(CHI(v) / v.size()) : T(0.);
}

template<typename T> inline Real MAG(const Vec<T> &v)
{
	return sqrt(real(DOT(v, v)));
}

template<typename T>
inline Vec<T> ABS(const Vec<T> &v)
{
	Vec<T> ret(v.size());
	for_Idx (i, 0, v.size()) ret[i] = ABS(v[i]);
	return ret;
}

VecReal ABS(const VecCmplx &v)
{
	VecReal ret(v.size());
	for_Idx(i, 0, v.size()) ret[i] = ABS(v[i]);
	return ret;
}

template<typename T>
inline Vec<T> EXP(const Vec<T> &v)
{
	Vec<T> ret(v.size());
	for_Idx(i, 0, v.size()) ret[i] = EXP(v[i]);
	return ret;
}

template<typename T>
inline Vec<T> SQR(const Vec<T> &a)
{
	Vec<T> ret = a;
	return ret *= ret;
}

template<typename T>
inline Vec<T> SQRT(const Vec<T> &a)
{
	Vec<T> ret(a.size());
	for_Idx(i, 0, a.size()) ret[i] = SQRT(a[i]);
	return ret;
}

template<typename T>
inline Vec<T> kronecker_product(const Vec<T> &a, const Vec<T> &b)
{
	Vec<T> ret(a.size() * b.size());
	for_Idx (i, 0, a.size()) for_Idx (j, 0, b.size()) ret[i * b.size() + j] = a[i] * b[j];
	return ret;
}

template<typename T>
inline VecReal real(const Vec<T> &v)
{
	VecReal r(v.size());
	for_Idx (i, 0, v.size()) r[i] = real(v[i]);
	return r;
}

template<typename T>
inline VecReal imag(const Vec<T> &v)
{
	VecReal r(v.size());
	for_Idx (i, 0, v.size()) r[i] = imag(v[i]);
	return r;
}

inline VecCmplx cmplx(const VecReal &a, const VecReal &b)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(a.size(), b.size());
#endif
	VecCmplx z(a.size());
	for_Idx (i, 0, z.size()) z[i] = cmplx(a[i], b[i]);
	return z;
}

inline VecCmplx cmplx(const VecReal &a)
{
	VecCmplx z(a.size());
	for_Idx (i, 0, z.size()) z[i] = cmplx(a[i], 0.);
	return z;
}

template<typename T>
inline Vec<T> INV(const Vec<T> &a)
{ 
	Vec<T> res(a.size());
	for_Int(i, 0, a.size()) res[i] = INV(a[i]);
	return res;
}

#endif /* _VEC_H_ */
