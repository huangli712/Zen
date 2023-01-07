#pragma once

/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

#include "densitymatrix.h"
#include "green.h"

// Krylov method for claculate the green function, for this methond is only need one space.

class Krylov : public DensityMat
{
private:
	const Real& ground_state_energy;
	const VecReal& ground_state_i;
	const Int position;
public:

private:
	// Using the Lanczos to mapping H to the tridiagonal matrix(step by step).
	void find_trdgnl_one_step(const VecReal& initial_vector, VecReal& v0, VecReal& v1, Real& a, Real& b, const SparseMatReal& sep_h);

public:
	Krylov(const DensityMat& old_one_emat_i, const Int position);

	ImGreen find_gf_from_krylov_space(const Int crtann);
};
