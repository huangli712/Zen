#pragma once

/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

#include "specs.h"
#include "prmtr.h"
#include "densitymatrix.h"
#include "nocspace.h"
#include "state.h"
#include "crrltfun.h"
#include "crrvec.h"
#include "krylov.h"
#include "impurity.h"

// impurity model

class NORG {
	// typedef VEC<VecInt> Tab;
	//*****************iteration***************
	Int iter_norg_cnt;						// count the NORG optimize iteration times.
	Int norg_stable_count;					// count the NORG stability iteration times.
	VecReal occnum_pre;						// occupation number record before one optimize iteration.
	VecReal occnum_lst;						// occupation number record after one optimize iteration.
	Real	groune_pre;						// ground state energy record before one optimize iteration.
	Real	groune_lst;						// ground state energy record after one optimize iteration.
	Real groundenererr;						// ground state energy Error: error between the two NORG iterations
	Real occupationerr;						// occupation Error			: error between the two NORG iterations
	Real correctionerr;						// correction Error			: error between two correction vector modify
	mutable VecCmplx checkg_pre;			// correction vec check point before one optimize iteration.
	mutable VecCmplx checkg_lst;			// correction vec check point after one optimize iteration.

public:
	const MyMpi &mm;						// parameters
	const Prmtr &p;							// parameters
	VEC<MatReal> uormat;					// unitary_orbital_rotation_matrix

	MatReal h0;
	VecReal final_ground_state;
	ImGreen impgreen;

	mutable NocSpace scsp;					// The main space.
	mutable DensityMat oneedm;				// The main space's density matrix.

	// mutable NocSpace scsp_1mone, scsp_1pone;
	// mutable Operator n1mone, n1pone;
	// const Tab& table_n1mone, table_n1pone;

private:
	// only change for first norg_set of the first div.
	VecInt nppso(const VecInt &a, Int positon)
	{
		VecInt nppso_i(a);
		Int i = abs(positon);
		if (positon > 0) nppso_i[(i - 1) * 2] += 1;
		if (positon < 0) nppso_i[(i - 1) * 2] -= 1;
		return nppso_i;
	}

	bool converged();

	bool green_modify_converged() const;
	bool green_modify_converged(Real correctionerr_i) const;

	void readmatrix(MatReal& m, const Str& file);

	// using the correct vector to modify the shortcut space.
	void upgrade_space(NocSpace& scsp_i, NocSpace& scsp_ipl, VecReal& state_pl, NocSpace& scsp_imi, VecReal& state_mi);

	VEC<MatReal> uormat_initialize();

	void show_the_nozero_number_of_tabel();
/*

	void update_by_correct_vec() const;

	void update_by_correct_vec(const Int orbital, const VecReal& ground_state_temp, const Real& ge,  const VEC<MatReal>& uormat_i) const;

	void rotation_the_space(const VEC<MatReal>& uormat_i) const;
*/
public:
	NORG(const MyMpi& mm_i, const Prmtr& prmtr_i);
	// NORG(const MyMpi& mm_i, const Impurity& imp_i, const Prmtr& prmtr_i);

	void up_date_h0_to_solve(const MatReal& h0_i);
	
/*
	// In this method we calculate the correction vector using the Krylov-space approach to modify the U.
	void modify_by_krylov_correct_vec();

	// (Deprecated) In this method we calculate the correction vector by the conjugategradient method.
	void modify_by_correct_vector() ;

	// (Deprecated)
	void writematrix(const MatReal& m, const Str U_name, Int iter_norg_cnt) const;
*/

//---------------------------------------calculate the physical operator---------------------------------

	Real sz_imp_sz_bath(const Int imp_postition, const VecReal& vgs_i);

	// To check the NO-interactions check.
	Real return_nointeractions_ground_state_energy(const MatReal& h0_i) const;
//--------------------------------------- for the Green function---------------------------------
	
	void get_gimp_by_krylov_CV_modify(Green& imp_i) const;
	// // Only use for test the validity.
	// void get_gimp_by_krylov_CV_modify(ReGreen& imp_i) const;

	Int get_gimp_with_possible_degeneracy(ImGreen& imp_i, Int iter_cont = 999);

	void get_gimp_by_krylov(const ImGreen& imp_i);

	void find_g(Green &g) ;
	
	Int find_g_with_possible_degeneracy(VEC<ReGreen> &g, Int kind = 0) ;

	// give the impurity green from krylov correction vec.
	void get_g_by_KCV(Green& imp_i);
	// give the impurity green matrix at onece from krylov correction vec.
	void get_g_by_KCV_spup(Green& imp_i);

	// git the impurity green from continue fraction.
	void get_g_by_CF(Green& imp_i);
	
//--------------------------------------- for the print out---------------------------------
	void write_norg_info(Int iter_cnt) const;

	void write_state_info(Int iter_cnt) const;
	
	// (Deprecated)
	MatReal save_transform_uormat();
};

