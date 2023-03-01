#pragma once
/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China)
date 2022-04-25
*/

#include "specs.h"
#include "prmtr.h"
#include "bath.h"
#include "model.h"
#include "green.h"
#include "impurity.h"
#include "norg.h"
#include "occler.h"



// A class for transform the data of ZEN.
// - Preprocessor. (APIzen)
//   - read hybrid function
//   - do the fitting
//   - Temporarily store the build matrix information required by NORG

class APIzen{
	const MyMpi& mm;					// parameters
	Prmtr& p;							// parameters
	Int num_orbital;					// The number of orbiatals.
	Int num_nondegenerate;				// The number of nondegenerate orbiatal.
	Real num_omg;						// Number of positive imaginary frequencies used for bath fitting.
	MatCmplx imfrq_hybrid_function;		// Imaginary frequencies hybrid function, nrows():number of nondegenerate orbiatal.
	VecInt or_deg_idx;					// orbitals  degenerate idx.
	VecReal solver_eimp_data;			// Impurity energy level for orbitals, correspondingly mean Impurity energy.
	// ctqmcdata solver_ctqmc_data;
	Real Uc, Jz, mu;
	Int nband;							
	Int norbs;

	
	// NORG coding console
	Int norm_mode, test_mode, fast_mode;
	VecInt restrain, distribute;

	// NORG test part
	//Int nimp;
	VecReal bathose, bathhop;

public:
	VEC<VecReal> t_ose;						// hopping integral for all sites
	VEC<VecReal> t_hyb;						// H_0 onset energy
	VecReal muvec;
	Int dmft_cnt;
private:
	void read_ZEN(const Str& file);

	// void fitting();

	void test_for_fitting(const Bath& bth, const ImGreen& hby_i, Int num = 666);
public:
	APIzen(const MyMpi& mm_i, Prmtr& p, const Str& file = empty_str, const Int mode = 0);

};