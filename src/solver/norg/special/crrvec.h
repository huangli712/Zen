#pragma once

/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

#include "operator.h"
#include "green.h"

// correct vector method is a way to get a vec \mid c(\omega+i\eta)\rangle

class Crrvec
{
	const MyMpi& mm;							// parameters
	const Prmtr& p;								// parameters
	const NocSpace& old_nosp;					// nature orbital space
	const NocSpace& new_nosp;					// nature orbital space
	const VecReal& ground_state;
	/*SparseMatReal minus_A;*/
	SparseMatReal sep_h;
	const Int crtorann;							// decide which green function will you calculate.
	Real ground_state_energy, omega, eta;

	const Int right_set, right_pos_in_div;

	//save the Krylov spcae method's a, b and .
	VecReal a_saved, b_saved, inner_m_saved;

public:
	const Operator& opr;						// The operator of the space.
	const VecReal ex_state;						// vector of ex-state in new nosp
	
	VecCmplx correct_vector;					// vector of correct vector

	MatCmplx correct_vecs;						// vector of correct vector for different check point

	VecCmplx correct_green;

private:
	// mapping the state by using the project operator.
	VecReal project_uplwer_parical_space(const VecReal &initial_vector, const Int crtann, const Int norg_set, const Int orbit_pos_in_div = 0) const;
	 
	// Int count_for_sign_in_one_div() const;

	SparseMatReal find_diagonal_sparse_matrix(Real number);

	SparseMatReal add_diagonal_sparse_matrix(VEC<VecInt> h_idx, Real ge0, Real omega);

	// Using the Lanczos to mapping H to the tridiagonal matrix(step by step).
	void find_trdgnl_one_step(const VecReal& initial_vector, VecReal& v0, VecReal& v1, Real& a, Real& b, const SparseMatReal& sep_h);

	// Using the tridiagonal Krylov space matirx to find each excition state.
	void find_excted_state_by_ab(const VecReal& initial_vector, VecReal& v0, VecReal& v1, VecReal& vec_a, VecReal& vec_b, Int& k, const SparseMatReal& sep_h);

	bool if_in_this_orbital(const VecOnb &exd_cf, const Int crtann, const Int norg_set, const Int orbit_pos_in_div) const;

public:
	//Crrvec(const MyMpi &mm_i, const Prmtr &prmtr_i, const NocSpace &old_nosp, const NocSpace &main_nosp, const VecReal &vgs_i, const Real &gs, const Int right_set, const Real omeg = 0., const Real eta_i = 0.01);
	Crrvec(const NocSpace &old_nosp, const Operator& main_opr, const VecReal &vgs_i, const Real &gs, const Int right_set_i, const Real omeg = 0., const Real eta_i = 0.01);

	Crrvec(const NocSpace &old_nosp_i, const Operator& main_opr, const VecReal &vgs_i, const Real &gs, const Int right_set_i, const VecReal &omega_point);
	Crrvec(const NocSpace &old_nosp_i, const Operator& main_opr, const VecReal &vgs_i, const Real &gs, const Int right_set_i, const VecCmplx &omega_point);
	Crrvec(const NocSpace &old_nosp_i, const Operator& main_opr, const Int right_set_i, const Int right_pos_in_div_i, const VecReal &vgs_i, const Real &gs);

	Crrvec(const NocSpace &old_nosp_i, Operator& main_opr, const VecReal &vgs_i, const Real &gs, Green& mat_green, const Int set_position = 0);

	VecReal operator*(const VecReal& KetVec);

	VecReal imag_to_real(const VecReal& KetVec);

	VecCmplx find_correct_vec();

	ImGreen find_gf();

	ImGreen find_gf_from_krylov_space();

	void prepare_Krylov_space();
	
	void krylov_update_state(VecCmplx omega_vec);

	VecCmplx krylov_space_for_green(Int left_set, Int left_pos_in_div,const VecCmplx& omega_vec);
	
	Vec<MatCmplx> krylov_space_for_green_matrix(const Green& mat_green);
};
