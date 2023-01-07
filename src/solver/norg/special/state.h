#pragma once

/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

#include "specs.h"
#include "prmtr.h"
#include "nocspace.h"	// cause now This class was only for the NORG imp solver, at the Shortcut space.

/// <summary>
/// In this class will find all the possible state(StateStatistics) for one gived Combined division (ComDiv).
/// impurity model
/// At the Shortcut space(NocSpace) we set first of two divisons as impurity and the active orbital(bath site).
/// For now it's support for find two fermion act on one state.
/// </summary>
class StateStatistics 
{
	typedef VEC<std::tuple<Int, Int, MatInt>> hopdata;
private:
	static const BnmlFast bf;
public:
	const NocSpace& space;								// parameters
	const Int div_idx;										// The divison idx at all Combined division space.
	const Int idx;											// The idx in one divison.
	const MatInt occ_n;										// occupy number in each division.
	// const MatInt sit_mat;								// spin - orbits number in each division.
	// const VecInt sit_n;									// spin - orbits number in each division.
	const ComDivs cfg;										// The configuration of the state.
	// const hopdata hopup;		// For the upper spin-orbits, a set of occupy number in each division
	// const hopdata hopdw;		// For the down spin-orbits, a set of occupy number in each division
	//ComDivs cfgdw;										// The configuration of the state.

	// const VecOnb cfgup;										// Here only for better reading reason.
	// const VecOnb cfgdw;										// Here only for better reading reason.

	// for the statistics of the state.
	// VEC<Int> filled_up;				// a set of been filled site number.
	// VEC<Int> filled_dw;				// a set of been filled site number.
	Vec< VEC<Int> > filled_spinless;
	// VEC<VecInt> off_diagonal_term;	// a set of valuable term with the annihilation, creation orbit's positon, and Colum idx(H_i).
private:
	ComDivs cfgdetail();
	// To statistic the occupation situation and filling the filled_spinless matrix.
	void statistics();



	// fermion anticommutativity, for two orbitals changed.
	//Int anticommu(const Int& pi,const Int& pj);

	// find the div's idx according to the "newdiv".
	Int find_newdiv_idx(MatInt newdiv);
public:
	//state(Int idxD, VecInt bases_i) :idx(idxD), bases(bases_i) {};
	StateStatistics(const Int& h_i, const Int& comdiv_idx, const NocSpace& s_i);
	// The out put<0>[i]: the annihilation divsion position; <1>[i]: the creation dision positon; <2>[i]: The new divs.
	hopdata divocchop_ingroup(const Int& ComDiv, Idx sets_n);
	// Output all off-diagonal H matrix term in one row.(For all sites hoping?)
	// [i][0]:annihilation orbit's position; [i][1]:creation orbit's positon; [i][2]:Colum idx(i); [i][3]:sign(fermion anticommutativity)
	// PS.:The matrix Max size: (N*N-(N-N_{imp})-N_{imp}*N_{imp},2), and it's spinless.
	VEC<VecInt> find_off_diagonal_term(const hopdata &hopup, const hopdata &hopdw);

	// Output all off-diagonal H matrix term with one spinless orbit.
	// [i][0]:annihilation orbit's position; [i][1]:creation orbit's positon; [i][2]:Colum idx(i); [i][3]:sign(fermion anticommutativity)
	VEC<VecInt> find_each_spiless_group_off_diagonal_term(const hopdata &hopspi, const Int group);
};
