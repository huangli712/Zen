/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2013 - 2017
*/

#ifndef _LLINTDOUBLE_H_
#define _LLINTDOUBLE_H_

#include "stdpfx.h"

const Int full_length_double_fraction = 53;
const Int half_length_double_fraction = 27;
const LLInt min_double_fraction = LLInt(1) << 52;
const LLInt max_double_fraction = (LLInt(1) << 53) - LLInt(1);

// represent double number d = f * 2 ^ e
// f is the fraction part of d and is a long long int number in [min_double_fraction, max_double_fraction] union {-1, 0, 1}
// the sign of d is incorporated into f
class LLIntDouble {
public:
	double d;
	LLInt f;
	Int e;
public:
	LLIntDouble() : d(0), f(0), e(0) {}
	LLIntDouble(double d_i): d(d_i)
	{
		if (ABS(d) == 0.) {
			f = 0;
			e = 0;
		} else if (ABS(d) >= max_double) {
			f = d < 0 ? -1 : 1;
			e = 2048;
		} else {
			double tmp = d;
			e = 0;
			while (ABS(tmp) > double(max_double_fraction)) { tmp /= 2; ++e; }
			while (ABS(tmp) < double(min_double_fraction)) { tmp *= 2; --e; }
			f = LLInt(tmp);
		}
	}
	Int check() const
	{
		double ff = double(f);
		Int ee = e;
		while (ee < 0) { ff /= 2.; ++ee; }
		while (ee > 0) { ff *= 2.; --ee; }
		if (ABS(ff) == 0. && ABS(d) == 0.)
			return 1;
		else if (ABS(ff) > max_double && ABS(d) > max_double)
			return 1;
		else if (ABS(ff) <= max_double && ABS(d) <= max_double && ABS(ff - d) / ABS(d) < 1.E-14)
			return 1;
		else {
			WRN(STR("LLIntDouble, d = ", d, ", f = ", f, ", e = ", e, ", ff = ", ff));
			return 0;
		}
	}
};

// f = full_length_double_fraction
// h = half_length_double_fraction
// (c2 << f) + c0 = a * b
// (binary length of x) <= f for x = a, b, c2, c0
// a, b, c2, c0 >= 0
inline void llintdouble_fraction_multiplication(const LLInt &a, const LLInt &b, LLInt &c2, LLInt &c0)
{
	const Int &half = half_length_double_fraction;
	const Int &full = full_length_double_fraction;
	if (a < 0 || b < 0 || a > max_double_fraction || b > max_double_fraction) ERR(STR("a = ", a, ", b = ", b));
	LLInt a0 = mod_by_mask(a, half);
	LLInt b0 = mod_by_mask(b, half);
	LLInt a1 = a >> half;
	LLInt b1 = b >> half;
	LLInt c1 = a1 * b0 + a0 * b1;
	LLInt c1_hgh = c1 >> half;
	LLInt c1_low = mod_by_mask(c1, half);
	c2 = a1 * b1 + c1_hgh;
	c0 = a0 * b0 + (c1_low << half);
	c2 = (c2 << (2 * half - full)) + (c0 >> full);
	c0 = mod_by_mask(c0, full);
}

#endif /* _LLINTDOUBLE_H_ */
