/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2013 - 2017
*/

#ifndef _TENSOR_H_
#define _TENSOR_H_

#include "stdpfx.h"
#include "vec.h"
#include "mat.h"

template<typename T>
// trns operates on rows of init, yields rows of dest
inline void tensor_transform(Mat<T> &dest, const Mat<T> &init, const Mat<T> &trns)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ3(dest.nrows(), init.nrows(), dest.ncols(), trns.nrows(), init.ncols(), trns.ncols());
#endif
	for_Idx (i, 0, init.nrows()) MUL(dest[i], trns, init[i]);
}

template<typename T>
// trns operates on the second dimension of tensor init, yields the second dimension of tensor dest
// init is of dimension s0 * n0 * s2
// dest is of dimension s0 * n1 * s2
inline void tensor_transform(Idx s0, Idx n1, Idx n0, Idx s2, 
							 Vec<T> &dest, const Vec<T> &init, const Mat<T> &trns)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ3_PS(n1, trns.nrows(), n0, trns.ncols(), s0 * n1 * s2 != dest.size(), s0 * n0 * s2 != init.size(), 
		NAVC8(s0, n1, n0, s2, dest.size(), init.size(), trns.nrows(), trns.ncols()));
#endif
	Vec<T> x(s2 * s0 * n0);
	Mat<T>(s2, s0 * n0, x.p()).tr(Mat<T>(s0 * n0, s2, init.p()));
	Vec<T> y(s2 * s0 * n1);
	tensor_transform(Mat<T>(s2 * s0, n1, y.p()), Mat<T>(s2 * s0, n0, x.p()), trns);
	Mat<T>(s0 * n1, s2, dest.p()).tr(Mat<T>(s2, s0 * n1, y.p()));
}

template<typename T>
// matrix on (om) operaters on dimension n (m) of init
// init is of dimension s0 * m0 * s1 * n0 * s2
// dest is of dimension s0 * m1 * s1 * n1 * s2
inline void tensor_transform(Idx s0, Idx m1, Idx m0, Idx s1, Idx n1, Idx n0, Idx s2, 
							 Vec<T> &dest, const Vec<T> &init, const Mat<T> &om, const Mat<T> &on)
{
	Vec<T> temp(s0 * m0 * s1 * n1 * s2);
	tensor_transform(s0 * m0 * s1, n1, n0, s2, temp, init, on);
	tensor_transform(s0, m1, m0, s1 * n1 * s2, dest, temp, om);
}

#endif /* _TENSOR_H_ */
