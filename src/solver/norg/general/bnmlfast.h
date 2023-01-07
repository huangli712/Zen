/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2012 - 2017
*/

#ifndef _BNMLFAST_H_
#define _BNMLFAST_H_

#include "stdpfx.h"
#include "vec.h"
#include "mat.h"
#include "number.h"

class BnmlFast {
private:
	static const Int N;		// condition: (M - 1) <= (N - 1) / 2
	static const Int M;		// 2147483647 = max_Int < 2333606220 = binomial(34, 17)
	static const Int S;
	static const Int L;		// binomial(65537, 2) = 2147516416 > 2147483647 = max_Int
	MatInt tm;				// of size M * N, binomial(0 : M - 1, 0 : N - 1), requiring (M - 1) <= (N - 1) / 2
	Vec<VecInt> tv;			// of size M * (n - N), binomial(S : M - 1, 0 : n - 1)
public:
	BnmlFast();
	Int operator()(Int n, Int m) const;		// (n < 0 || m < 0 || m > n) && binomial(n, m) < max_Int is prerequisite
	~BnmlFast() {}
};

#endif /* _BNMLFAST_H_ */
