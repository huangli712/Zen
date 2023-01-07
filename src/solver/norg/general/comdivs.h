#pragma once

/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2021-2022
*/

#include "stdpfx.h"
#include "bnmlfast.h"
#include "onb.h"
#include "vec.h"

//
// A state formed under a specific set of Combined division(ComDivs) and with the idx for this combined division.
// Hilbert subspace (in one one gived Combined division(ComDivs) ) 
// labeled by (I,C^N1_M1,C^N2_M2,C^N3_M3,...C^Nn_Mn) = (I, vec N, vec M) = (idx, VecInt ne, VecInt ns)


class ComDivs {
private:
	static const BnmlFast bf;
public:
	typedef std::tuple<VEC<VecInt>, VEC<VecInt>> OrbTuple;
	const VecInt ne;				// a set of occupy number in each division 
	const VecInt ns;				// a set of spin-orbits number in each division 
	Int idx;						// The idx for one Combined division
	const VecOnb cf;				// a set of configuration. 
	VEC<VecInt> div_orb_e;			// each division contain the electron-orbital's Idx.
	VEC<VecInt> div_orb_c;			// each division contain the cavity-orbital's Idx.
	Int FreeBaseCode(VecInt& rep, VecInt& base) const;
	VecInt FreeBaseDeCode(Int idx, VecInt& base) const;
	Int count_electrons(const Int& lw, const Int& up) const;
public:
	// A quick way to bilid the ComDivs, and only support for find the cf.
	ComDivs(Int i_i, VecInt ne_i, VecInt ns_i, bool quick);

	// A quick way to bilid the ComDivs, and only support for find the cf.
	ComDivs(Int i_i, MatInt ne_i, MatInt ns_i, bool quick);

	// according the idx construct the ComDivs, get the VecOnb cf.
	ComDivs(Int i_i, VecInt ne_i, VecInt ns_i);

	// according the idx construct the ComDivs, get the VecOnb cf.
	ComDivs(Int i_i, MatInt ne_i, MatInt ns_i);

	// according the newcf construct the ComDivs, get the idx.
	ComDivs(VecOnb& newcf, VecInt ne_i, VecInt ns_i);
	
	// according the newcf construct the ComDivs, get the idx.
	ComDivs(VecOnb& newcf, MatInt ne_i, MatInt ns_i);

	// Only suit for fermion sign changed by c_i^+ c_j or c_j^+ c_i
	Int sgn(Int pi, Int pj) const {
		if (pi > pj) SWAP(pi, pj);
		return 1 - 2 * (1 & count_electrons(pi + 1, pj));
	}

	// ----Function---

	// return subspace index for a specific ns
	VecOnb cf_ComDivs(const Int idx_in);
	// return the div_orb with two index, [i]:the i div number; [i][j]: in i div,the j-th orbit, and the value is the orbit idx.
	OrbTuple find_div_orb();
	Int find_cf_idx(const VecOnb& newcf);
	// 

};




