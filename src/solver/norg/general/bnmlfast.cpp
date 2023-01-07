/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2012 - 2017
*/

#include "bnmlfast.h"

const Int BnmlFast::N = 34;			// condition: (M - 1) <= (N - 1) / 2
const Int BnmlFast::M = 17;			// 2147483647 = max_Int < 2333606220 = binomial(34, 17)
const Int BnmlFast::S = 3;
const Int BnmlFast::L = 65537;		// binomial(65537, 2) = 2147516416 > 2147483647 = max_Int

BnmlFast::BnmlFast()
{
	tm.reset(M, N);
	for_Int (i, 0, M) {
		for (Int j = 0; j < i; ++j) tm[i][j] = 0;
		for (Int j = i; j < N; ++j) tm[i][j] = binomial(j, i, 0);
	}

	tv.reset(M);
	for_Int (i, S, M) {
		Int n = 0;
		for_Int (j, N, L) {
			if (binomial(j, i, 0) < 0) {
				n = j;
				break;
			}
		}
		tv[i].reset(n - N);
		for_Int (j, N, n) tv[i][j - N] = binomial(j, i, 0);
	}
}

Int BnmlFast::operator()(Int n, Int m) const
{
#ifdef _ASSERTION_
	if (n < 0 || m < 0 || m > n) ERR(STR("BnmlFast()(", n, ", ", m, ")"));
#endif
	if (m > (n >> 1)) m = n - m;
	if (n < N) return tm[m][n];
	if (m == 1) return n;
	if (m == 0) return 1;
#ifdef _ASSERTION_
	if (m == 2 && n >= L) return -1;
#endif
	if (m == 2) return (n & 1) ? n * ((n - 1) >> 1) : (n >> 1) * (n - 1);
#ifdef _ASSERTION_
	if (m >= M) return -1;
	if (n >= tv[m].size() + N) return -1;
#endif
	return tv[m][n - N];	
}
