/*
code developed by
	Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2020-10-26
*/

#ifndef _BIJECTION_H_
#define _BIJECTION_H_

#include "stdpfx.h"
#include "vec.h"

// class for bijection: Idx -> T
// from 0, 1, ..., n - 1 to objects of type T
// bijection[i = 0, 1, .., n - 1] = object of type T

template<typename T>
class Bijection {
private:
	Idx n;		// number of elements
	Vec<T> idx2t;
	std::map<T, Idx> t2idx;
public:
	Bijection() : n(0) {}
	void init(const std::vector<T>& v) {
		n = v.size();
		idx2t = v;
		for_Idx(i, 0, n) {
			t2idx[v[i]] = i;
		}
		if (t2idx.size() < n) {
			ERR(STR("there are duplicates in v"));
		}
	}
	inline const T& obj(Idx i) const {
		return idx2t[i];
	}
	inline Idx idx(const T& t) const {
		return t2idx.at(t);
	}
	inline bool exist(const T& t) const {
		if (t2idx.find(t) == t2idx.end()) {
			return false;
		}
		else {
			return true;
		}
	}
	Idx size() const { return n; }
};

#endif /* _BIJECTION_H_ */
