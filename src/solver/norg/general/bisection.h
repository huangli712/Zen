/*
code developed by
        Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2015 - 2017
*/

#ifndef _BISECTION_H_
#define _BISECTION_H_

#include "stdpfx.h"
#include "vec.h"
#include "mat.h"

template<typename T>
// v is in ascending order
// if (v.size() <= 1) return 0;
// else if (key < v[0]) return 0;
// else if (key > v[v.size() - 1]) return v.size() - 1;
// else return l; such that v[l] == key || v[l] < key && key < v[l + 1]
inline Idx bisection(const Vec<T> &v, T key)
{
	Idx l = 0;
	Idx r = v.size();
	while (l + 1 < r) {
		Idx m = (l + r) >> 1;			// l < m < r
		if (key > v[m]) l = m;
		else if (key < v[m]) r = m;
		else return m;
	}
	return l;
}

#endif /* _BISECTION_H_ */
