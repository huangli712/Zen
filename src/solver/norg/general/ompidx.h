/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2013 - 2017
*/

#ifndef _OMPIDX_H_
#define _OMPIDX_H_

// occupancy number basis

#include "stdpfx.h"
#include "vec.h"
#include "slctsort.h"

// occupancy number basis state (for example |01001101>)
// a binary number, 1 for occupied, 0 for unoccupied, spinless
// 0-based numbering is used
// the state of the i-th site is denoted by the i-th digit (from low to high)
// Int pos denote the position of an site
// for ns <= 32, use UInt for bs
// for ns <= 64, use ULLInt for bs
// for ns > 64, change the code
// (const Int &lw, const Int &up) means [lw, up)

class OmpIdx {
private:
	Idx n;					// number of indices
	VEC<Int> idx;			// index
	VEC<Int> cst;			// cost
public:
	OmpIdx(): n(0) {}
	void push_back(const Int &i, const Int &c) {
		idx.push_back(i);
		cst.push_back(c);
		++n;
	}
	void descend() { slctsort(cst, idx, std::greater<Int>()); }
	Int operator[](const Idx &ii) const { return idx[ii]; }
	Int i(const Idx &ii) const { return idx[ii]; }
	Vec<Int> idxvec() { return Vec<Int>(idx); }
	Idx size() const { return n; }
};

class OmpIdxIdx {
private:
	Idx n;					// number of indices
	VEC<Int> idx0;			// index
	VEC<Int> idx1;			// index
	VEC<Int> cst;			// cost
public:
	OmpIdxIdx(): n(0) {}
	void push_back(const Int &i0, const Int &i1, const Int &c) {
		idx0.push_back(i0);
		idx1.push_back(i1);
		cst.push_back(c);
		++n;
	}
	void descend() { slctsort(cst, idx0, idx1, std::greater<Int>()); }
	std::pair<Int, Int> operator[](const Idx &ii) const { return std::make_pair(idx0[ii], idx1[ii]); }
	Int i0(const Idx &ii) const { return idx0[ii]; }
	Int i1(const Idx &ii) const { return idx1[ii]; }
	Vec<Int> idx0vec() { return Vec<Int>(idx0); }
	Vec<Int> idx1vec() { return Vec<Int>(idx1); }
	Idx size() const { return n; }
};

#endif /* _OMPIDX_H_ */
