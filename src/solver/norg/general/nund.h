#pragma once

/*
coded by Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2021-02-26
*/

#include "stdpfx.h"
#include "bnmlfast.h"

// Hilbert subspace labeled by (Nu, Nd) = (nu, nd) = (ne_up, ne_dw)
class NuNd {
private:
	static const BnmlFast bf;
	Int ns;		// number of lattice sites
	Int nd;		// number of spin-down electrons
	Int nu;		// number of spin-up electrons
public:
	NuNd(Int ns_i, Int nd_i, Int nu_i) : ns(ns_i), nd(nd_i), nu(nu_i) {}
	NuNd(Int ns_i, Idx idx_nund) : ns(ns_i), nd(idx_nund % (ns + 1)), nu(idx_nund / (ns + 1)) {}
	bool operator<(const NuNd& b) const {
		if (ns == b.ns && nu == b.nu) return nd < b.nd;
		return ns == b.ns ? nu < b.nu : ns < b.ns;
	}
	// return subspace index for a specific ns
	Idx idx_nund() const { return nd + nu * (ns + 1); }
	Idx dim_dw() const { return bf(ns, nd); }
	Idx dim_up() const { return bf(ns, nu); }
	Idx dim() const { return dim_dw() * dim_up(); }
	bool exist() const {
		return 0 <= nd && nd <= ns && 0 <= nu && nu <= ns && 0 <= ns;
	}
	NuNd nd_add_one() const {
		return NuNd(ns, nd + 1, nu);
	}
	NuNd nd_sub_one() const {
		return NuNd(ns, nd - 1, nu);
	}
	NuNd nu_add_one() const {
		return NuNd(ns, nd, nu + 1);
	}
	NuNd nu_sub_one() const {
		return NuNd(ns, nd, nu - 1);
	}
	Int get_ne_dw() const { return nd; }
	Int get_ne_up() const { return nu; }
	Int get_ns() const { return ns; }
	// idx_basis = idx_dw + idx_up * dim_dw()
	Idx idx_dw(Idx idx_basis) const {
		return idx_basis % dim_dw();
	}
	// idx_basis = idx_dw + idx_up * dim_dw()
	Idx idx_up(Idx idx_basis) const {
		return idx_basis / dim_dw();
	}
	// idx_basis = idx_dw + idx_up * dim_dw()
	Idx idx_basis(Idx idx_dw, Idx idx_up) const {
		return idx_dw + idx_up * dim_dw();
	}
};
