#pragma once

/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

#include "operator.h"
#include "green.h"

// Correlation function The method is easily generalized to calculate correlations 
// $$\left\langle A^{\dagger}(t) B\left(t^{\prime}\right)\right\rangle$$ for any pair of observables $$A^{+}$$and $$B$$ 
// and, once the ground-state properties $$\left(\left|\psi_{0}\right\rangle\right,and\left.E_{0}\right)$$ are known,
// all correlations can be obtained with a small computational effort.
// Here \psi_{0} mean ground state(gs), ground state energy(e_g).
// DOI:: Physical review letters 59.26 (1987): 2999.

class CrrltFun : public Operator
{
	typedef std::tuple<VecReal, VecReal> Vectrdgnl;
	typedef std::tuple<StdVecReal, StdVecReal> VECtrdgnl;
	const NocSpace& old_nosp;				// nature orbital space
	const NocSpace& new_nosp;				// nature orbital space
	//VecReal vgs;							// vector ground state(gs)
	const Int crtorann;
	//const Real& gse;					// ground state energy(gse)	
	//const Int ex_orbital_position;

	//save the continue fraction a and b.
	VecReal a_saved, b_saved;

public:
	const VecReal ex_state;				// vector of ex-state in new scsp

private:
	//  (Deactivate) Using the Lanczos to find the tridiagonal matrix.
	Vectrdgnl find_trdgnl(const VecReal& initial_vector, const Int crtann);


	// Using the Lanczos to mapping H to the tridiagonal matrix(First step).
	//VECtrdgnl find_trdgnl_first(const VecReal& initial_vector);

	// Using the Lanczos to mapping H to the tridiagonal matrix(step by step).
	void find_trdgnl_one_step(const VecReal& initial_vector, VecReal& v0, VecReal& v1, Real& a, Real& b, const SparseMatReal& sep_h);

	// mapping the state by using the project operator.
	VecReal project_uplwer_parical_space(const VecReal &initial_vector, const Int crtann, const Int norg_set, const Int orbit_pos_in_div = 0) const;

	bool if_in_this_orbital(const VecOnb &exd_cf, const Int crtann, const Int norg_set, const Int orbit_pos_in_div) const;

public:
	CrrltFun(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& old_nosp_i, const NocSpace& main_nosp, const VecReal& vgs_i, const Int position = 0);
	CrrltFun(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& old_nosp_i, const NocSpace& main_nosp, const Tab& table, const VecReal& vgs_i, const Int position = 0);
	CrrltFun(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& old_nosp_i, const NocSpace& main_nosp, const Tab& table, const VecReal& vgs_i, const Int pos_in_set, const Int pos_in_div);

	// CrrltFun(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& old_nosp_i, const NocSpace& main_scsp, const Tab& table, const MatReal& vdegs_i, const Int position = 0);
	// ? CrrltFun(const NocSpace& old_nosp_i, const Operator& main_op, const VecReal& vgs_i, const Int position = 0);

	// // (Deactivate) Caluate the density-density correlation function.
	// $$G_{A}(Z)=\left\langle\psi_{0}\left|A^{\dagger}(Z-H)^{-1}A\right|\psi_{0}\right\rangle,(Z=\omega+i\eta+E_{0})$$
	ImGreen find_density_density_correlation_function(const Real &ge0);

	void find_gf_greater(const Real& ge0, Green &g0);

	void find_gf_lesser(const Real& ge0, Green &g0);

	
	void find_gf_greater(const Real& ge0, Green &g0, Int kind);

	void find_gf_lesser(const Real& ge0, Green &g0, Int kind);
};
