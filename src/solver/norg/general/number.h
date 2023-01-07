/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2012 - 2017
*/

#ifndef _NUMBER_H_
#define _NUMBER_H_

#include "stdpfx.h"

template<typename T>
Int count_binary_ones(T x)
{
	Int sum = 0;
	while (x) {
		++sum;
		x &= x - 1;
	}
	return sum;
}

template<typename T>
inline T GCD(T a, T b)
{
	if (a < 0 || b < 0 || a == 0 && b == 0) ERR(STR("GCD(", a, ", ", b, ")"));
	if (a == 0) return b;
	if (b == 0) return a;
	Int d = 0;
	while (!(a & 1) && !(b & 1)) {
		a >>= 1;
		b >>= 1;
		++d;
	}
	while (a != b) {
		if (!(a & 1)) a >>= 1;
		else if (!(b & 1)) b >>= 1;
		else if (a > b) a = (a - b) >> 1;
		else b = (b - a) >> 1;
	}
	return a << d;
}

template<typename T>
inline T LDM(T a, T b)
{
	if (a < 0 || b < 0 || a == 0 && b == 0) ERR(STR("LDM(", a, ", ", b, ")"));
	return a / GCD(a, b) * b;
}

template<typename T>
inline Int binary_trailing_zeros(T a)
{
	ULLInt b = ULLInt(a);
	if (!b) return 64;
	Int i = 0;
	while (!(b & 1)) {
		b >>= 1;
		++i;
	}
	return i;
}

template<typename T>
inline Int binary_trailing_ones(T a)
{
	ULLInt b = ULLInt(a);
	Int i = 0;
	while (b & 1) {
		b >>= 1;
		++i;
	}
	return i;
}

// 20! < max_LLInt = 9223372036854775807 < 21!
// n < 0 || n > 20 is not allowed
inline LLInt factorial(Int n, Int wrn = 1)
{
	if (n < 0) ERR(STR("factorial(", n, ")"));
	if (n > 20) {
		WRN(STR("factorial(", n, ") > max_LLInt = ", max_LLInt));
		return -1;
	}
	LLInt a = 1;
	for (Int i = 1; i <= n; ++i) {
		a *= i;
	}
	return a;
}

// return -1 if binomial(n, m) > max_int = 2147483647 < 2333606220 = binomial(34, 17)
// (n < 0 || m < 0 || m > n) is not allowed
inline Int binomial(Int n, Int m, Int wrn = 1)
{
#ifdef _ASSERTION_
	if (n < 0 || m < 0 || m > n) ERR(STR("binomial(", n, ", ", m, ")"));
#endif
	if (m > n / 2) m = n - m;
	if (m > 20) {
		if (wrn) WRN(STR("binomial(", n, ", ", m, ")", " > (137846528820 = binomial(40, 20))"));
		return -1;
	}
	if (m == 0) return 1;
	if (m == 1) return n;

	// now 4 <= n && 2 <= m && m <= n / 2

	const Int max_Int = 2147483647;
	const LLInt max_LLInt = 9223372036854775807;
	const LLInt overflow = max_LLInt / n;
	Int l = (n + n - m + 1) / 2;

	LLInt a = 1;
	for (Int i = n - m + 1; i <= l; ++i) {
		if (a > overflow) {
			if (wrn) WRN(STR("binomial(", n, ", ", m, "), i = ", i, ", overflow = ", overflow));
			return -1;
		}
		a *= i;
	}

	LLInt b = 1;
	for (Int i = l + 1; i <= n; ++i) {
		if (b > overflow) {
			if (wrn) WRN(STR("binomial(", n, ", ", m, "), i = ", i, ", overflow = ", overflow));
			return -1;
		}
		b *= i;
	}

	LLInt d = 1;
	for (Int i = 1; i <= m; ++i) {
		d *= i;
	}

	LLInt g;
	g = GCD(a, d);
	a /= g;
	d /= g;

	g = GCD(b, d);
	b /= g;
	d /= g;
	
#ifdef _ASSERTION_
	if (d != 1) ERR(STR("binomial(", n, ", ", m, "), d = ", d, ", a = ", a, ", b = ", b));
#endif

	if (a > max_LLInt / b) {
		if (wrn) WRN(STR("binomial(", n, ", ", m, ") > max_LLInt = ", max_LLInt, ", a = ", a, ", b = ", b));
		return -1;
	}

	a *= b;

	if (a > max_Int) {
		if (wrn) WRN(STR("binomial(", n, ", ", m, ") = ", a, " > max_Int = ", max_Int));
		return -1;
	}

	return Int(a);
}

#endif /* _NUMBER_H_ */
