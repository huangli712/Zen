/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2013 - 2017
*/

#ifndef _HPSI_H_
#define _HPSI_H_

#include "stdpfx.h"
#include "vec.h"
#include "mat.h"

// class for sparse hamiltonian
// row(i, VecIdx &rc, Vec<Th> &rh) should guarantee the uniqueness of element indices in a row

/* -----------------------
   NOT TESTED AND NOT USED
   ----------------------- */

template<typename Th, typename Tv, typename F>
class HPsi {
private:
	Idx n;				// dimension
	Vec<Vec<Th> > h;	// sparse matrix
	Vec<VecIdx> c;		// column index
public:
	typedef Th Type_h;									// make Th available externally
	typedef Th Type_v;									// make Tv available externally
	HPsi(Idx n_i, F &row);								// row(i, VecIdx &rc, Vec<Th> &rh) generates row i in Vec r
	void operator()(const Vec<Tv> &x, Vec<Tv> &y);		// y = H x
	void hermitianize();								// hermitianize
	Int ishermitian();									// if h is hermitian
	Idx element_address(Idx i, Idx j);					// return k such that c[i][k] == j;
};

template<typename Th, typename Tv, typename F>
HPsi<Th, Tv, F>::HPsi(Idx n_i, F &row): n(n_i), h(n), c(n)
{
	for_Idx (i, 0, n) {
		VecIdx rc;
		Vec<Th> rh;
		row(i, rc, rh);
#ifdef _CHECK_DIMENSION_MATCH_
		if (rc.size() != rh.size()) ERR(NAV4(i, n, rc.size(), rh.size()));
#endif
		sort_Vec_Vec(rc, rh, std::less<Idx>());
		c[i].reset(rc.size());
		h[i].reset(rh.size())
		c[i] = rc;
		h[i] = rh;
	}
}

template<typename Th, typename Tv, typename F>
void HPsi<Th, Tv, F>::operator()(const Vec<Tv> &x, Vec<Tv> &y)
{
#ifdef _CHECK_DIMENSION_MATCH_
	if (n != x.size() || n != y.size()) ERR(NAV3(n, x.size(), y.size()));
#endif
	for_Idx (i, 0, n) {
		Tv tmp = Tv(0);
		for_Idx (j, 0, h[i].size()) tmp += h[i][j] * x[c[j]];
		y[i] = tmp;
	}
}

template<typename Th, typename Tv, typename F>
Int HPsi<Th, Tv, F>::ishermitian()
{
	Real max_abs_element = 0.;
	for_Idx (i, 0, n) for_Idx (j, 0, h[i].size()) {
			if (max_abs_element < ABS(h[i][j])) max_abs_element = ABS(h[i][j]);
	}
	if (max_abs_element < sqrt(eps_Real)) max_abs_element = sqrt(eps_Real);
	Real eps = 100 * sqrt(Real(n)) * eps_Real * max_abs_element;

	for_Idx (row, 0, n) for_Idx (j, 0, h[i].size()) {
		Idx col = c[i][j];
		Idx k = element_address(col, row);
		if (k >= n) {
			WRN("not hermitian, element not found" + NAVC6(n, row, col, j, k, h[row][j]));
			return 0;
		}
		if (ABS(h[row][j] - cnjg(h[col][k])) > eps) {
			WRN("not hermitian" + NAVC9(n, row, col, j, k, c[col][k], h[row][j], h[col][k], eps));
			return 0;
		}
	}
	return 1;
}

template<typename Th, typename Tv, typename F>
void HPsi<Th, Tv, F>::hermitianize()
{
	for_Idx (row, 0, n) for_Idx (j, 0, h[i].size()) {
		Idx col = c[i][j];
		Idx k = element_address(col, row);
		if (k >= n) {
			h[row][j] = 0.;
		} else {
			Th tmp = 0.5 * (h[row][j] + cnjg(h[col][k]));
			h[row][j] = tmp;
			h[col][k] = cnjg(tmp);
		}
	}
}

template<typename Th, typename Tv, typename F>
inline Idx HPsi<Th, Tv, F>::element_address(Idx i, Idx j)
{
#ifdef _CHECK_BOUNDS_
	if (i >= n || j >= n) ERR(NAV3(i, j, n));
#endif
	for_Idx (k, 0, c[i].size()) if (c[i][k] == j) return k;
	return n;
}

#endif /* _HPSI_H_ */
