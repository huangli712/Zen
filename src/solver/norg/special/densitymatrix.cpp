/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

#include "densitymatrix.h"
using namespace std;

DensityMat::DensityMat(const MyMpi& mm_i, const Prmtr& prmtr_i, NocSpace& scsp_i) :Operator(mm_i, prmtr_i, scsp_i)
{
}
MatReal DensityMat::one_electron_density_matrix(Int wish_nev)
{
	MatReal multi_oedm(p.norbit, p.norbit, 0.);		// one - electron density matrix at the Multi - states.
	// if (mm) WRN("lowest_eigpairs BEGIN");
	// if (mm)WRN("find_one_electron_density_matrix BEGIN"+NAV(ground_states.size()));
	multi_oedm = find_one_electron_density_matrix(lowest_eigpairs(scsp.dim, true,  wish_nev));
	return multi_oedm;
}

VEC<MatReal> DensityMat::find_one_electron_density_matrix(const MatReal& state, const Tab& table_i) 
{
	// if (mm) WRN("find_one_electron_density_matrix BEGIN: ");
	VecReal state_eff(state.ncols(), 0.);
	for_Int(i, 0, state.nrows()) state_eff += state[i];
	state_eff.normalize();
	VecPartition row_H(mm.np(), mm.id(), scsp.dim);

	VEC<MatReal> D_splited;
	for_Int(i, 0, p.nI2B.size()) {
		MatReal temp(p.nI2B[i] + scsp.sit_mat[i][0], p.nI2B[i] + scsp.sit_mat[i][0], 0.);
		D_splited.push_back(std::move(temp));
	}

	const Tab& hop_op(table_i);
	for_Idx(pos, 0, hop_op[2].size())
	{
		if (abs(hop_op[2][pos]) <= scsp.hopint.size())
		{
			Int row_n((abs(hop_op[2][pos]) - 1) / p.norbit), col_n((abs(hop_op[2][pos]) - 1) % p.norbit); // "1" mean the real_length = abs(hop_op[2][pos]) - 1, since to distinguish the minus sign.
			Int check_point(0);
			for_Int(i, 0, p.nI2B.size()) {
				check_point = 0;
				if (row_n >= 0 && row_n < (p.nI2B[i] + scsp.sit_mat[i][0]) && col_n >= 0 && col_n < (p.nI2B[i] + scsp.sit_mat[i][0]))
				{
					Int coefficient_comm = hop_op[2][pos] >= 0 ? 1 : -1;
					Real amplitude = coefficient_comm * state_eff[hop_op[0][pos] + row_H.bgn()] * state_eff[hop_op[1][pos]];
					D_splited[i][row_n][col_n] += amplitude;
					check_point++;
					break;
				}
				else {
					row_n -= p.nI2B[i] + scsp.sit_mat[i][0];
					col_n -= p.nI2B[i] + scsp.sit_mat[i][0];
				}
			}
#ifdef _ASSERTION_
			if (check_point == 0)ERR("Not in the hopint but you save this bases(Many body base), why? For D mapping");
#endif
		}
	}
	VEC<MatReal> D;
	for_Int(i, 0, p.norg_sets) {
		MatReal D_i(mm.Allreduce(D_splited[i]));
		D.push_back(std::move(D_i));
	}
#ifdef _ASSERTION_
	for (const auto& d : D)
		if (d.if_symmetric() == false)ERR("The one_electron_density_matrix is not hermitian!!");
#endif
	return std::move(D);
}

VEC<MatReal> DensityMat::correct_one_electron_density_matrix(const VecReal& state, Crrvec& corstate_p, Crrvec& corstate_m, const VecReal& omega_point)
{
	VecReal state_eff = state;
	state_eff.normalize();
	VecPartition row_H(mm.np(), mm.id(), scsp.dim);

	VEC<MatReal> D_splited;
	for_Int(i, 0, p.nI2B.size()) {
		MatReal temp(p.nI2B[i] + scsp.sit_mat[i][0], p.nI2B[i] + scsp.sit_mat[i][0], 0.);
		D_splited.push_back(std::move(temp));
	}
	// for the ground state's D.
	const Tab& hop_op(table);
	for_Idx(pos, 0, hop_op[2].size())
	{
		if (abs(hop_op[2][pos]) <= scsp.hopint.size())
		{
			Int row_n((abs(hop_op[2][pos]) - 1) / p.norbit), col_n((abs(hop_op[2][pos]) - 1) % p.norbit); // "1" mean the real_length = abs(hop_op[2][pos]) - 1, since to distinguish the minus sign.
			Int check_point(0);
			for_Int(i, 0, p.nI2B.size()) {
				check_point = 0;
				if (row_n >= 0 && row_n < (p.nI2B[i] + scsp.sit_mat[i][0]) && col_n >= 0 && col_n < (p.nI2B[i] + scsp.sit_mat[i][0]))
				{
					Int coefficient_comm = hop_op[2][pos] >= 0 ? 1 : -1;
					Real amplitude = coefficient_comm * state_eff[hop_op[0][pos] + row_H.bgn()] * state_eff[hop_op[1][pos]];
					D_splited[i][row_n][col_n] += amplitude;
					check_point++;
					break;
				}
				else {
					row_n -= p.nI2B[i] + scsp.sit_mat[i][0];
					col_n -= p.nI2B[i] + scsp.sit_mat[i][0];
				}
			}
#ifdef _ASSERTION_
			if (check_point == 0)ERR("Not in the hopint but you save this bases(Many body base), why? For D mapping");
#endif
		}
	}

	// for_Int(j, 0, omega_point.size()){
		corstate_p.krylov_update_state( I * cmplx(omega_point));
		// for the corstate_p's D.
		find_density_matrix_by_Crrvec(D_splited, corstate_p);
		corstate_m.krylov_update_state( I * cmplx(omega_point));
		// for the corstate_m's D.
		find_density_matrix_by_Crrvec(D_splited, corstate_m);
	// }

	VEC<MatReal> D;
	for_Int(i, 0, p.norg_sets) {
		MatReal D_i(mm.Allreduce(D_splited[i]));
		D.push_back(std::move(D_i));
	}
#ifdef _ASSERTION_
	for (const auto& d : D)
		if (d.if_symmetric() == false)ERR("The one_electron_density_matrix is not hermitian!!");
#endif
	return std::move(D);
}


VEC<MatReal> DensityMat::correct_one_electron_density_matrix(const VecReal& state, Crrvec& corstate1_p, Crrvec& corstate1_m, Crrvec& corstate2_p, Crrvec& corstate2_m)
{
	VecReal state_eff = state;
	state_eff.normalize();
	VecPartition row_H(mm.np(), mm.id(), scsp.dim);

	VEC<MatReal> D_splited;
	for_Int(i, 0, p.nI2B.size()) {
		MatReal temp(p.nI2B[i] + scsp.sit_mat[i][0], p.nI2B[i] + scsp.sit_mat[i][0], 0.);
		D_splited.push_back(std::move(temp));
	}
	// for the ground state's D.
	const Tab& hop_op(table);
	for_Idx(pos, 0, hop_op[2].size())
	{
		if (abs(hop_op[2][pos]) <= scsp.hopint.size())
		{
			Int row_n((abs(hop_op[2][pos]) - 1) / p.norbit), col_n((abs(hop_op[2][pos]) - 1) % p.norbit); // "1" mean the real_length = abs(hop_op[2][pos]) - 1, since to distinguish the minus sign.
			Int check_point(0);
			for_Int(i, 0, p.nI2B.size()) {
				check_point = 0;
				if (row_n >= 0 && row_n < (p.nI2B[i] + scsp.sit_mat[i][0]) && col_n >= 0 && col_n < (p.nI2B[i] + scsp.sit_mat[i][0]))
				{
					Int coefficient_comm = hop_op[2][pos] >= 0 ? 1 : -1;
					Real amplitude = coefficient_comm * state_eff[hop_op[0][pos] + row_H.bgn()] * state_eff[hop_op[1][pos]];
					D_splited[i][row_n][col_n] += amplitude * 0.5;
					check_point++;
					break;
				}
				else {
					row_n -= p.nI2B[i] + scsp.sit_mat[i][0];
					col_n -= p.nI2B[i] + scsp.sit_mat[i][0];
				}
			}
#ifdef _ASSERTION_
			if (check_point == 0)ERR("Not in the hopint but you save this bases(Many body base), why? For D mapping");
#endif
		}
	}

	
	// for the corstate_p's D.
	find_density_matrix_by_Crrvec(D_splited, corstate1_p);
	find_density_matrix_by_Crrvec(D_splited, corstate2_p);

	// for the corstate_m's D.	
	find_density_matrix_by_Crrvec(D_splited, corstate1_m);
	find_density_matrix_by_Crrvec(D_splited, corstate2_m);

	VEC<MatReal> D;
	for_Int(i, 0, p.norg_sets) {
		MatReal D_i(mm.Allreduce(D_splited[i]));
		D.push_back(std::move(D_i));
	}
#ifdef _ASSERTION_
	for (const auto& d : D)
		if (d.if_symmetric() == false)ERR("The one_electron_density_matrix is not hermitian!!");
#endif
	return std::move(D);
}


VEC<MatReal> DensityMat::correct_one_electron_density_matrix(const VecReal& state, const Crrvec& corstate_p, const Crrvec& corstate_m)
{
	VecReal state_eff = state;
	state_eff.normalize();
	VecPartition row_H(mm.np(), mm.id(), scsp.dim);
	// Real d_weight_x    = 1.6; 	//4times
	Real d_weight_x = 1.; 		//2times
	// Real d_weight_x = (4./7.); 	//1times

	VEC<MatReal> D_splited;
	for_Int(i, 0, p.nI2B.size()) {
		MatReal temp(p.nI2B[i] + scsp.sit_mat[i][0], p.nI2B[i] + scsp.sit_mat[i][0], 0.);
		D_splited.push_back(std::move(temp));
	}
	// for the ground state's D.
	const Tab& hop_op(table);
	for_Idx(pos, 0, hop_op[2].size())
	{
		if (abs(hop_op[2][pos]) <= scsp.hopint.size())
		{
			Int row_n((abs(hop_op[2][pos]) - 1) / p.norbit), col_n((abs(hop_op[2][pos]) - 1) % p.norbit); // "1" mean the real_length = abs(hop_op[2][pos]) - 1, since to distinguish the minus sign.
			Int check_point(0);
			for_Int(i, 0, p.nI2B.size()) {
				check_point = 0;
				if (row_n >= 0 && row_n < (p.nI2B[i] + scsp.sit_mat[i][0]) && col_n >= 0 && col_n < (p.nI2B[i] + scsp.sit_mat[i][0]))
				{
					Int coefficient_comm = hop_op[2][pos] >= 0 ? 1 : -1;
					Real amplitude = coefficient_comm * state_eff[hop_op[0][pos] + row_H.bgn()] * state_eff[hop_op[1][pos]];
					D_splited[i][row_n][col_n] += amplitude * 0.25 *d_weight_x;
					check_point++;
					break;
				}
				else {
					row_n -= p.nI2B[i] + scsp.sit_mat[i][0];
					col_n -= p.nI2B[i] + scsp.sit_mat[i][0];
				}
			}
#ifdef _ASSERTION_
			if (check_point == 0)ERR("Not in the hopint but you save this bases(Many body base), why? For D mapping");
#endif
		}
	}

	// for the corstate_p's D.
	find_density_matrix_by_Crrvec(D_splited, corstate_p);

	// for the corstate_m's D.
	find_density_matrix_by_Crrvec(D_splited, corstate_m);

	VEC<MatReal> D;
	for_Int(i, 0, p.norg_sets) {
		MatReal D_i(mm.Allreduce(D_splited[i]));
		D.push_back(std::move(D_i));
	}
#ifdef _ASSERTION_
	for (const auto& d : D)
		if (d.if_symmetric() == false)ERR("The one_electron_density_matrix is not hermitian!!");
#endif
	return std::move(D);
}

void DensityMat::find_density_matrix_by_Crrvec(VEC < MatReal>& D_splited, const Crrvec& corstate_i)
{
	VecPartition row_H(mm.np(), mm.id(), corstate_i.opr.scsp.dim);
	MatReal corrvec_real = real(corstate_i.correct_vecs); for_Int(i, 0, corrvec_real.nrows()) corrvec_real[i].normalize();
	MatReal corrvec_imag = imag(corstate_i.correct_vecs); for_Int(i, 0, corrvec_imag.nrows()) corrvec_imag[i].normalize();
	VecReal oneparticalex = corstate_i.ex_state; oneparticalex.normalize();
	// Real d_weight_y    = 0.8; 	//4times
	Real d_weight_y = 1.; 		//2times
	// Real d_weight_y = (8./7.); 	//1times

	const Tab& hop_op(corstate_i.opr.table);
	for_Idx(pos, 0, hop_op[2].size())
	{
		if (abs(hop_op[2][pos]) <= corstate_i.opr.scsp.hopint.size())
		{
			Int row_n((abs(hop_op[2][pos]) - 1) / p.norbit), col_n((abs(hop_op[2][pos]) - 1) % p.norbit); // "1" mean the real_length = abs(hop_op[2][pos]) - 1, since to distinguish the minus sign.
			Int check_point(0);
			for_Int(i, 0, p.nI2B.size()) {
				check_point = 0;
				if (row_n >= 0 && row_n < (p.nI2B[i] + scsp.sit_mat[i][0]) && col_n >= 0 && col_n < (p.nI2B[i] + scsp.sit_mat[i][0]))
				{
					Int coefficient_comm = hop_op[2][pos] >= 0 ? 1 : -1;
					for_Int(j, 0, corstate_i.correct_vecs.nrows()){
						Real amplitude = coefficient_comm * corrvec_real[j][hop_op[0][pos] + row_H.bgn()] * corrvec_real[j][hop_op[1][pos]];
							amplitude += coefficient_comm * corrvec_imag[j][hop_op[0][pos] + row_H.bgn()] * corrvec_imag[j][hop_op[1][pos]];
							amplitude += coefficient_comm * oneparticalex[hop_op[0][pos] + row_H.bgn()] * oneparticalex[hop_op[1][pos]];
						D_splited[i][row_n][col_n] += d_weight_y * 0.25 * amplitude * INV(corstate_i.correct_vecs.nrows() * 2.);// fist 2: G^gatter, G^lesser; second 2: Rell and imag part.
					}
					check_point++;
					break;
				}
				else {
					row_n -= p.nI2B[i] + scsp.sit_mat[i][0];
					col_n -= p.nI2B[i] + scsp.sit_mat[i][0];
				}
			}
#ifdef _ASSERTION_
			if (check_point == 0)ERR("Not in the hopint but you save this bases(Many body base), why? For D mapping");
#endif
		}
	}
}

Real DensityMat::sum_off_diagonal() const{
	Real sum = 0 ;
	auto test_dm(dm);
	for_Int(i, 0, p.norg_sets) if(i%2 == 0) {
		MatReal test_dm_bath = test_dm[i].truncate(scsp.sit_mat[i][0], scsp.sit_mat[i][0], p.nI2B[i] + scsp.sit_mat[i][0], p.nI2B[i] + scsp.sit_mat[i][0]);
		for_Int(r, 0, test_dm_bath.nrows()) test_dm_bath[r][r] -= test_dm_bath[r][r];
		VecReal vectemp(test_dm_bath.vec());
		// if(mm) WRN(NAV(test_dm_bath));
		sum += SQRT(vectemp.norm_sqr() / (vectemp.size() - test_dm_bath.nrows()));
	}
	return sum/(p.norg_sets / 2);
}

VEC<MatReal> DensityMat::find_unitary_orbital_rotation_matrix()
{
	VEC<MatReal> bathdm;
	for_Int(i, 0, p.norg_sets) {
		bathdm.push_back(dm[i].truncate(1, 1, p.nI2B[i] + 1, p.nI2B[i] + 1));
		if(mm) WRN(NAV(dm[i][0][0]))
	}
	// if (mm) WRN(NAV3(dm[0], dm[1], dm[2]));

	VEC<VecReal> evalue;
	for_Int(i, 0, p.norg_sets) {
		VecReal evalu_i(bathdm[i].nrows(), 0.);
		// if(mm) WRN("New dm for bath" + NAV3(i, bathdm[i], bathdm.size()));
		heevr(bathdm[i], evalu_i); bathdm[i] = bathdm[i].tr(); evalue.push_back(std::move(evalu_i));
		// if(mm) WRN("New U for bath" + NAV2(bathdm[i], evalue[i]));
		for_Int(j, 0, (p.nI2B[i] / 2.)) {
			SWAP(bathdm[i][j], bathdm[i][p.nI2B[i] - j - 1]);
			SWAP(evalue[i][j], evalue[i][p.nI2B[i] - j - 1]);
		}
		//DBG("New uorm111" + NAV3(i, bathdm[i], evalue[i]));
	}
	if(mm) WRN(NAV3(evalue[0].mat(1,p.nI2B[0]), evalue[1].mat(1,p.nI2B[1]), evalue[2].mat(1,p.nI2B[2])));
	for_Int(i, 0, bathdm.size()) bathdm[i] = bathdm[i - (i%2)];

	occupationnumber = evalue;
	return bathdm;
}

VEC<VecReal> DensityMat::check_dm_get_occupation_number() const{
	auto test_dm(dm);
	VEC<VecReal> occupation_nubmer;
	for_Int(i, 0, p.norg_sets) {
		VecReal evalu_i(test_dm[i].nrows(), 0.);
		heevr(test_dm[i], evalu_i); test_dm[i] = test_dm[i].tr();
		if(mm) WRN(NAV3(i,evalu_i,SUM(evalu_i)));
		if (i % 2 == 0) occupation_nubmer.push_back(evalu_i);
	}
	return occupation_nubmer;
}

// (Deactivate) In Here the state a row is a state, and the number of row means degeneracy
MatReal DensityMat::find_one_electron_density_matrix(const MatReal& state)
{
//#ifdef _CHECK_DIMENSION_MATCH_
//	ASSERT_EQ(state.ncols(), degeneracy);
//#endif
	// if (mm)WRN("find_one_electron_density_matrix BEGIN: ");
	VecReal state_eff(state.ncols(), 0.);
	for_Int(i, 0, state.nrows()) state_eff += state[i];
	state_eff.normalize();
	// if (mm)WRN("find_one_electron_density_matrix BEGIN: " + NAV2(state_eff.isnormalized(), SUM(state_eff)) );
	VecPartition row_H(mm.np(), mm.id(), scsp.dim);
	MatReal D_splited(p.norbit, p.norbit, 0.);
	for_Int(h_i, row_H.bgn(), row_H.end())
	{
		StateStatistics a(h_i, scsp.wherein_NocSpace(h_i), scsp);
		// For diagonally terms.
		for (const auto &x : a.filled_spinless)
		{
			for (const auto &i : x)
			{
				D_splited[i][i] += state_eff[h_i] * state_eff[h_i];
			}
		}
		vector<VecInt> off_dt;
		for_Int(idx_sets, 0, p.norg_sets){
			vector<VecInt> off_dt_next(a.find_each_spiless_group_off_diagonal_term(a.divocchop_ingroup(a.div_idx, idx_sets),idx_sets));
			off_dt. insert(off_dt. end(), off_dt_next. begin(), off_dt_next. end());
		}
		for_Int(i, 0, off_dt.size()) {
			D_splited[off_dt[i][1]][off_dt[i][0]] += off_dt[i][3] * state_eff[h_i] * state_eff[off_dt[i][2]];
		}
	}
	MatReal D(mm.Allreduce(D_splited));
	if (mm)WRN("find_one_electron_density_matrix FINISHED!!! " + NAV(D));
#ifdef _ASSERTION_
	if (D.if_symmetric() == false)ERR("The one_electron_density_matrix is not hermitian!!");
#endif
// #ifdef _ASSERTIONFORVSCODE_
// 	if (D.if_symmetric() == false)ERR("The one_electron_density_matrix is not hermitian!!");
// #endif
	return D;
}

// (Deactivate) for the OrthonormalityRecover.
MatReal DensityMat::OrthonormalityRecover(const MatReal& mat_i)
{
	MatReal mat(mat_i.tr());
	for_Int(i, 0, mat.nrows()) {
		for_Int(j, 0, i) {
			mat[i] -= mat[j] * DOT(mat[i], mat[j]);
		}
		mat[i].normalize();
	}
	return mat.tr();
}
