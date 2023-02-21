#pragma once

/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

#include "operator.h"
#include "crrvec.h"

//impurity model

class DensityMat : public Operator
{
public:
	VEC<MatReal> dm;				// density matrix
	// VEC<MatCmplx> dm_cmplx;		// density matrix
	// occupation number: for each row stand the bath's occupation for one specific imp orbit.
	VEC<VecReal> occupationnumber;

private:
	//(Deactivate) To build the one-electron density matrix at the Multi-states.
	MatReal one_electron_density_matrix(Int wish_nev = 1);

	// for all the orbital in one matrix.(Deactivate)
	// find the one-electron density matrix at ground state in shortcut space.
	// and the one-electron density matrix is splited storage in different prosses.
	MatReal find_one_electron_density_matrix(const MatReal& state);

	// find the one-electron density matrix by table.
	VEC<MatReal> find_one_electron_density_matrix(const MatReal& state, const Tab& table);

	// find the one-electron density matrix by table within one effective state.
	VEC<MatReal> correct_one_electron_density_matrix(const VecReal& state, const Crrvec& corstate_p, const Crrvec& corstate_m);

	// find the one-electron density matrix by table within the excitation state's krylov space.
	VEC<MatReal> correct_one_electron_density_matrix(const VecReal& state, Crrvec& corstate_p, Crrvec& corstate_m, const VecReal& omega_point);

	// find the one-electron density matrix by table within the excitation state's krylov space.
	VEC<MatReal> correct_one_electron_density_matrix(const VecReal& state, Crrvec& corstate1_p, Crrvec& corstate1_m, Crrvec& corstate2_p, Crrvec& corstate2_m);

	// find the one-electron density matrix by correct vector.
	void find_density_matrix_by_Crrvec(VEC < MatReal>& D_splited, const Crrvec& corstate_i);

public:
	DensityMat(const MyMpi& mm_i, const Prmtr& prmtr_i, NocSpace& scsp_i);
	DensityMat(const MyMpi& mm_i, const Prmtr& prmtr_i, NocSpace& scsp_i, Str tab_name);

	// Only diagonalize the bath orbital, but reture the whole orbital unitary orbital rotation matrix U
	VEC<MatReal> find_unitary_orbital_rotation_matrix();

	// To update calculate density matrix.
	void update() { dm = find_one_electron_density_matrix(lowest_eigpairs(scsp.dim), table); }

	// To update density matrix by corrstate.
	void update(const Crrvec& corstate_p, const Crrvec& corstate_m) { dm = correct_one_electron_density_matrix(ground_state, corstate_p, corstate_m); }

	// To update calculate density matrix by input state.
	void update(const VecReal& state, const Crrvec& corstate_p, const Crrvec& corstate_m) {
		dm = correct_one_electron_density_matrix(state, corstate_p, corstate_m);
	}

	// To update calculate density matrix by input Krylov sapce, with omega vec (tow orbital version).
	void update(const VecReal& state, Crrvec& corstate1_p, Crrvec& corstate1_m, Crrvec& corstate2_p, Crrvec& corstate2_m) {
		// state = lowest_eigpairs(scsp.dim)[0];
		dm = correct_one_electron_density_matrix(state, corstate1_p, corstate1_m, corstate2_p, corstate2_m);
	}

	// // (Deprecated) To update calculate density matrix by input Krylov sapce, with omega vec.
	// void update(const VecReal& state, Crrvec& corstate_p, Crrvec& corstate_m, const VecReal& omega_point) {
	// 	dm = correct_one_electron_density_matrix(state, corstate_p, corstate_m, omega_point);
	// }
	
	//
	Real sum_off_diagonal() const;
	// (Deactivate) for the OrthonormalityRecover.
	MatReal OrthonormalityRecover(const MatReal &mat);

	VEC<VecReal> check_dm_get_occupation_number() const;
};
