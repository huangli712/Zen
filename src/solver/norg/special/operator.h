#pragma once

/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

#include "specs.h"
#include "prmtr.h"
#include "nocspace.h"
#include "state.h"

// impurity model
typedef Vec<VEC<Int>> Tab;

class Operator {
public:
	const MyMpi& mm;				// parameters
	const Prmtr& p;					// parameters
	const NocSpace& scsp;			// NocSpace

	Tab table;						// The movement by two femi operator contain.
	Real groundstate_energy;		// The ground state energy on this shortcut restratin.

	
	VecReal ground_state;			// if degeneracy happened, we only get one of them.
private:

public:
	// Expand the Shortcut space under Number of particles(NumberSpa).
	Operator(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& s_i);
	Operator(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& s_i, const Tab &per_table);
	Operator(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& s_i, Str tab_name);
	
	// Operator::Operator(const Operator &per):
	// mm(per.mm), p(per.p), scsp(per.scsp), table(per.table) , 
	// groundstate_energy(per.groundstate_energy), ground_state(per.ground_state){}
	
	// [i][0]row Idx; [i][1]colum Idx;[i][2] for the position for the hopint.
	Tab find_h_idx();

	void clear(){ for_Int(i, 0, 3) VEC<int> ().swap(table[i]);}

	// (Deactivate) Find the Hamiltonian operator at shortcut space, separated storage in different processes.
	SparseMatReal find_hmlt();
	// By using the find_h_idx() to find the idx first then speed up the procedure.
	SparseMatReal find_hmlt(const Tab h_idx) const;


	// Using Lanczos algorithm to find several lowest eigenpairs of a Hermitian matrix
	MatReal lowest_eigpairs(const Idx n, bool if_need_fast = true, Int wish_nev = 1);

	//Compute the Single particle excited state(sn_prtcl_ex_state) 
	//here only support for the div[0] 
	//$$\left|f_{0}\right\rangle=A\left|\psi_{0}\right\rangle$$
	VecReal sn_prtcl_ex_state(const Int ex_orbital_position, const VecReal ground_state, const Int crtann)const;

	//(deprecated)
	VecReal particle_number_Inner_product(const Int position, const Int crtann)const;	

	void save_the_Tab(Tab& tab, Str name) const;
	Tab read_the_Tab(Str name) const;
};