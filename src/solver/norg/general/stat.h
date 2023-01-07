/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2017-08
*/

#ifndef _STAT_H_
#define _STAT_H_

#include "stdpfx.h"
#include "vec.h"

#define statf_Real(fin, fex) inline Real fin(const Real &x) { return fex(x); }
#define statf_Cmplx(fin, fex) inline Cmplx fin(const Cmplx &x) { return cmplx(fex(real(x)), fex(imag(x))); }
#define statf_VecReal(fin, fex) inline VecReal fin(const VecReal &x) { return fex(x); }
#define statf_VecCmplx(fin, fex) inline VecCmplx fin(const VecCmplx &x) { return cmplx(fex(real(x)), fex(imag(x))); }

namespace StatF
{
	statf_Real(sqr, SQR)
	statf_Real(abs, ABS)
	statf_Real(sqrt, SQRT)

	statf_Cmplx(sqr, SQR)
	statf_Cmplx(abs, ABS)
	statf_Cmplx(sqrt, SQRT)

	statf_VecReal(sqr, SQR)
	statf_VecReal(abs, ABS)
	statf_VecReal(sqrt, SQRT)

	statf_VecCmplx(sqr, SQR)
	statf_VecCmplx(abs, ABS)
	statf_VecCmplx(sqrt, SQRT)
}

#undef statf_Real
#undef statf_Cmplx
#undef statf_VecReal
#undef statf_VecCmplx

// avg + dev statistics, avg = average, dev = deviation
template<typename T>
class StatVec {
private:
	Int count;
	Vec<T> sum;
	Vec<T> sum_sqr;
public:
	StatVec(Idx size = 0) : count(0), sum(size, 0.), sum_sqr(size, 0.) {}
	void reset(Idx size) {
		count = 0;
		sum.reset(size, 0.);
		sum_sqr.reset(size, 0.);
	}
	void accumulate(const Vec<T> &entry) {
		++count;
		sum += entry;
		sum_sqr += StatF::sqr(entry);
	}
	Vec<T> avg() {
		return sum * T(INV(count));
	}
	Vec<T> var() {
		const Vec<T>& chi_sqr = sum_sqr * T(INV(count)) - StatF::sqr(avg());
		return StatF::abs(chi_sqr) * T(Real(count) / Real(count - 1));
	}
	Vec<T> dev() {
		return StatF::sqrt(var());
	}
	Vec<T> dev_of_avg() {
		return dev() * T(SQRT(INV(count)));
	}
};

// avg + dev statistics, avg = average, dev = deviation
template<typename T>
class Stat {
private:
	Int count;
	T sum;
	T sum_sqr;
public:
	Stat() : count(0), sum(0.), sum_sqr(0.) {}
	void accumulate(const T &entry) {
		++count;
		sum += entry;
		sum_sqr += StatF::sqr(entry);
	}
	T avg() {
		return sum * T(INV(count));
	}
	T var() {
		const T& chi_sqr = sum_sqr * T(INV(count)) - StatF::sqr(avg());
		return StatF::abs(chi_sqr) * T(Real(count) / Real(count - 1));
	}
	T dev() {
		return StatF::sqrt(var());
	}
	T dev_of_avg() {
		return dev() * T(SQRT(INV(count)));
	}
};

typedef StatVec<Real> StatVecReal;
typedef StatVec<Cmplx> StatVecCmplx;

#endif /* _STAT_H_ */
