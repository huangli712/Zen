/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2017-08
*/

#ifndef _STTSTCS_H_
#define _STTSTCS_H_

#include "stdpfx.h"
#include "vec.h"

#define sttstcsf_Real(fin, fex) inline Real fin(const Real &x) { return fex(x); }
#define sttstcsf_Cmplx(fin, fex) inline Cmplx fin(const Cmplx &x) { return cmplx(fex(real(x)), fex(imag(x))); }
#define sttstcsf_VecReal(fin, fex) inline VecReal fin(const VecReal &x) { return fex(x); }
#define sttstcsf_VecCmplx(fin, fex) inline VecCmplx fin(const VecCmplx &x) { return cmplx(fex(real(x)), fex(imag(x))); }

namespace SttstcsF
{
	sttstcsf_Real(sqr, SQR)
	sttstcsf_Real(abs, ABS)
	sttstcsf_Real(sqrt, SQRT)

	sttstcsf_Cmplx(sqr, SQR)
	sttstcsf_Cmplx(abs, ABS)
	sttstcsf_Cmplx(sqrt, SQRT)

	sttstcsf_VecReal(sqr, SQR)
	sttstcsf_VecReal(abs, ABS)
	sttstcsf_VecReal(sqrt, SQRT)

	sttstcsf_VecCmplx(sqr, SQR)
	sttstcsf_VecCmplx(abs, ABS)
	sttstcsf_VecCmplx(sqrt, SQRT)
}

// avg + dev statistics, avg = average, dev = deviation
template<typename T>
class SttstcsVec {
private:
	Int count;
	Vec<T> sum;
	Vec<T> sum_sqr;
public:
	SttstcsVec(Idx size = 0) : count(0), sum(size, 0.), sum_sqr(size, 0.) {}
	void reset(Idx size) {
		count = 0;
		sum.reset(size, 0.);
		sum_sqr.reset(size, 0.);
	}
	void accumulate(const Vec<T> &entry) {
		++count;
		sum += entry;
		sum_sqr += SttstcsF::sqr(entry);
	}
	Vec<T> avg() {
		return sum * T(INV(count));
	}
	Vec<T> dev() {
		Vec<T> chisqr = sum_sqr * T(INV(count)) - SttstcsF::sqr(avg());
		return SttstcsF::sqrt(SttstcsF::abs(chisqr));
	}
};

// avg + dev statistics, avg = average, dev = deviation
template<typename T>
class Sttstcs {
private:
	Int n;
	T sum;
	T sum_sqr;
public:
	Sttstcs() : n(0), sum(0.), sum_sqr(0.) {}
	void accumulate(const T &entry) {
		++n;
		sum += entry;
		sum_sqr += SttstcsF::sqr(entry);
	}
	T avg() {
		return sum * INV(n);
	}
	T dev() {
		T chisqr = sum_sqr * INV(n) - SttstcsF::sqr(avg());
		return SttstcsF::sqrt(SttstcsF::abs(chisqr));
	}
};

typedef SttstcsVec<Real> SttstcsVecReal;
typedef SttstcsVec<Cmplx> SttstcsVecCmplx;

#endif /* _STTSTCS_H_ */
