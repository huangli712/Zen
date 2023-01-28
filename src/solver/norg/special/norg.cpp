/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

#include "norg.h"
#define condition iter_norg_cnt < p.iter_max_norg  && !converged()

NORG::NORG(const MyMpi& mm_i, const Prmtr& prmtr_i) :
	mm(mm_i), p(prmtr_i), impgreen(prmtr_i.norbs, prmtr_i), uormat(uormat_initialize()),
	occupation_err(1.), energy_err(1.), correctionerr(1.), h0(prmtr_i.norbit, prmtr_i.norbit, 0.), occnum_pre(SUM(prmtr_i.nI2B), 0.),
	iter_norg_cnt(0), groune_pre(1.e99), groune_lst(0.), occnum_lst(SUM(prmtr_i.nI2B), 0.),
	scsp(prmtr_i, h0, prmtr_i.npartical), oneedm(mm, prmtr_i, scsp),norg_stable_count(0)
{
	show_the_nozero_number_of_tabel();
	// WRN(NAV(STR("tste")));
	if(mm) WRN(NAV(scsp.dim));
}

NORG::NORG(const MyMpi& mm_i, const Prmtr& prmtr_i, VecInt nparticals) :
	mm(mm_i), p(prmtr_i), impgreen(prmtr_i.norbs, prmtr_i), uormat(uormat_initialize()),
	occupation_err(1.), energy_err(1.), correctionerr(1.), h0(prmtr_i.norbit, prmtr_i.norbit, 0.), occnum_pre(SUM(prmtr_i.nI2B), 0.),
	iter_norg_cnt(0), groune_pre(1.e99), groune_lst(0.), occnum_lst(SUM(prmtr_i.nI2B), 0.),
	scsp(prmtr_i, h0, nparticals), oneedm(mm, prmtr_i, scsp),norg_stable_count(0)
{
	show_the_nozero_number_of_tabel();
	// WRN(NAV(STR("tste")));
	if(mm) WRN(NAV(scsp.dim));
}

// NORG::NORG(const MyMpi& mm_i, const Impurity& imp_i, const Prmtr& prmtr_i) :
// 	mm(mm_i), p(prmtr_i), impgreen(prmtr_i.norbs, prmtr_i), uormat(uormat_initialize()),
// 	occupation_err(1.), energy_err(1.), correctionerr(1.), h0(imp_i.h0), occnum_pre(SUM(prmtr_i.nI2B), 0.),
// 	iter_norg_cnt(0), groune_pre(1.e99), groune_lst(0.), occnum_lst(SUM(prmtr_i.nI2B), 0.),
// 	scsp(prmtr_i, imp_i.h0, prmtr_i.npartical), oneedm(mm, prmtr_i, scsp), norg_stable_count(0)
// {
// 	// up_date_h0_to_solve(h0);
// 	// get_gimp(impgreen);
// }

void NORG::up_date_h0_to_solve(const MatReal& h0_i) {
	if(mm) std::cout << std::endl;						// blank line
	h0 = h0_i;
	// if (mm) PIO(NAV2(h0,scsp.dim));
	scsp.coefficient = scsp.set_row_primeter_by_gived_mat(uormat, h0_i);
	oneedm.update(); final_ground_state = oneedm.ground_state;
	// if(mm) PIO(NAV(oneedm.sum_off_diagonal()));
	groune_lst = oneedm.groundstate_energy;
	if (mm)	{
		scsp.print();
		// write_norg_info(iter_norg_cnt);
		// write_state_info(iter_norg_cnt);
	}
	while (iter_norg_cnt < p.iter_max_norg && !converged()) {
		iter_norg_cnt++;
		if(mm) PIO("The iteration countting: " + NAV(iter_norg_cnt));
		VEC<MatReal> uormat_new(oneedm.find_unitary_orbital_rotation_matrix());
		for_Int(i, 0, uormat.size()) uormat[i] = uormat_new[i] * uormat[i];
		// if(mm) WRN(NAV(uormat[0]));
		scsp.coefficient = scsp.set_row_primeter_by_gived_mat(uormat, h0_i);	//if (mm) scsp.print();
		// if (mm)PIO("ground_state size" + NAV(oneedm.ground_state.size()));
		groune_pre = groune_lst;	occnum_pre = occnum_lst;
		oneedm.update(); final_ground_state = oneedm.ground_state;
		// if(mm) PIO(NAV(oneedm.sum_off_diagonal()));
		occnum_lst = VECVectoVec(oneedm.occupationnumber);
		groune_lst = oneedm.groundstate_energy;
		// if(mm) WRN(NAV(oneedm.dm[0]));
		if (mm) {
			// WRN(NAV(iter_norg_cnt));
			// for_Int(i, 0, uormat.size()) WRN(NAV(uormat_new[i]));
			// write_norg_info(iter_norg_cnt);
			// write_state_info(iter_norg_cnt);
		}
		occupation_err = SQRT(SUM(SQR(occnum_pre - occnum_lst)) / p.norbit);
		energy_err = 2 * (groune_pre - groune_lst) / (ABS(groune_pre) + ABS(groune_lst) + 1.);
		if (mm) {
			PIO(NAV2(energy_err, occupation_err));
			std::cout << "groundE_pre " << iofmt("sci") << groune_pre << std::setw(4) << "   groundE_lst " << groune_lst  << "  " << present() << std::endl;
		}
		if(mm) std::cout << std::endl;						// blank line
	}
	iter_norg_cnt = 0;		occupation_err = 1.;		energy_err = 1.;
	final_ground_state = oneedm.ground_state;	norg_stable_count = 0.;
	
	// Finish the NORGiteration.
}

void NORG::get_g_by_KCV(Green& imp_i)
{
	NocSpace scsp_1pone(p, h0, nppso(p.npartical, 1));
	NocSpace scsp_1mone(p, h0, nppso(p.npartical, -1));
	Operator n1pone(mm, p, scsp_1pone);
	Operator n1mone(mm, p, scsp_1mone);
	for_Int(i, 0, imp_i.g[0].nrows()) {
	// for_Int(i, 0, 3){
		// {Int i(0);
			scsp_1pone.coefficient = scsp_1pone.set_row_primeter_by_gived_mat(uormat, h0);
			scsp_1mone.coefficient = scsp_1mone.set_row_primeter_by_gived_mat(uormat, h0);

			Crrvec greaer(scsp, n1pone, 0, i, final_ground_state, groune_lst);
			Crrvec lesser(scsp, n1mone, 0, i, final_ground_state, groune_lst);
			
			{// test
				if(mm) WRN(NAV2(greaer.ex_state.norm(), lesser.ex_state.norm()));
			}

		for_Int(j, 0, imp_i.g[0].ncols()) {
		// for_Int(j, 0, 3){
			// {Int j(i);
			VecCmplx green_function, vec_z;
			if(imp_i.type_info() == STR("ImGreen")) vec_z = p.Im_z;
			else if(imp_i.type_info() == STR("ReGreen")) vec_z = p.Re_z;
			else ERR("The imp_i.type_info() was wrong");
			green_function  = greaer.krylov_space_for_green(0, j, vec_z);
			green_function += lesser.krylov_space_for_green(0, j, vec_z);
			for_Int(n, 0, imp_i.g.size()) imp_i[n][i][j] = green_function[n];
		}
	if(mm) std::cout << present() << std::endl;
	}
}

void NORG::get_g_by_KCV_spup(Green& imp_i)
{
	StdVecInt difference = {1, -1};
	for(const auto ii: difference)
	{
		NocSpace scsp_1(p, h0, nppso(p.npartical, ii));
		Operator opr_sub(mm, p, scsp_1);
		scsp_1.coefficient = scsp_1.set_row_primeter_by_gived_mat(uormat, h0);
		Crrvec greaer(scsp, opr_sub, final_ground_state, groune_lst, imp_i);
	}
	// {
	// 	NocSpace scsp_1mone(p, h0, nppso(p.npartical, -1));
	// 	Operator n1mone(mm, p, scsp_1mone);
	// 	scsp_1mone.coefficient = scsp_1mone.set_row_primeter_by_gived_mat(uormat, h0);
	// 	Crrvec lesser(scsp, n1mone, final_ground_state, groune_lst, imp_i);
	// }

	if(mm) std::cout << present() << std::endl;
}

void NORG::get_gimp(Green& imp_i)
{
	for_Int(i, 0, p.nband) {
		StdVecInt difference = {(i+1), -(i+1)};
		for(const auto ii: difference)
		{
			NocSpace scsp_sub(p, h0, nppso(p.npartical, ii));
			Operator opr_sub(mm, p, scsp_sub);
			scsp_sub.coefficient = scsp_sub.set_row_primeter_by_gived_mat(uormat, h0);
			CrrltFun temp_green(mm, p, scsp, scsp_sub, opr_sub.table, final_ground_state, i * 2);
			ImGreen green_function(1, p);
			if(ii > 0) temp_green.find_gf_greater(groune_lst, green_function);
			if(ii < 0) temp_green.find_gf_lesser(groune_lst, green_function);
			// green_function.write("continued_fraction-imp_green_function", 2 * i + 1);
			for_Int(n, 0, green_function.nomgs) imp_i[n][i][i] += green_function[n][0][0];
		}
		if (mm) PIO("finished the " + STR(i) + " find_g_norg   " + present());
	}
}

void NORG::readmatrix(MatReal& m, const Str& file)
{
	IFS ifs(file);
	if (!ifs) {
		if (mm) WRN(STR("file opening failed with ") + NAV(file))
	}
	else {
	for_Int(i, 0, m.nrows()) {
		for_Int(j, 0, m.ncols()) ifs >> m[i][j];
	}
	}
	ifs.close();
}

void NORG::get_g_by_CF(Green& imp_i)
{
	NocSpace scsp_1pone(p, h0, nppso(p.npartical, 1));
	NocSpace scsp_1mone(p, h0, nppso(p.npartical, -1));
	Operator n1pone(mm, p, scsp_1pone);
	Operator n1mone(mm, p, scsp_1mone);

	{//n=0 using the lanczos continue fraction.
		scsp_1pone.coefficient = scsp_1pone.set_row_primeter_by_gived_mat(uormat, h0);
		scsp_1mone.coefficient = scsp_1mone.set_row_primeter_by_gived_mat(uormat, h0);

		for_Int(i, 0, imp_i.g[0].nrows()){
			CrrltFun greaer(mm, p, scsp, scsp_1pone, n1pone.table, final_ground_state, 0 * 2, i);
			CrrltFun lesser(mm, p, scsp, scsp_1mone, n1mone.table, final_ground_state, 0 * 2, i);

			{// test
				if(mm) WRN(NAV2(greaer.ex_state.norm(), lesser.ex_state.norm()));
			}
			ImGreen green_function(1, p);
			greaer.find_gf_greater(groune_lst, green_function);
			lesser.find_gf_lesser(groune_lst, green_function);
			// green_function.write("continued_fraction-imp_green_function", 2 * i + 1);
			for_Int(n, 0, green_function.nomgs) imp_i[n][i][i] = green_function[n][0][0];
			if (mm) PIO("finished the find_g_norg");

		}
	}
}

VEC<MatReal> NORG::uormat_initialize()
{
	VEC<MatReal> uormat_i;
	for_Int(i, 0, p.nI2B.size()) {
		MatReal temp(dmat(p.nI2B[i], 1.));
		uormat_i.push_back(std::move(temp));
	}
	return std::move(uormat_i);
}

bool NORG::converged() 
{
		const Real gseerr = ABS(energy_err);
		if (gseerr > 1.E-3) { return false; }
		else if (gseerr < 1.E-10) { return true; }
		else if (energy_err < 1.E-5) { 
			norg_stable_count++;
			if (norg_stable_count > 10) return true;
			}
		else if (iter_norg_cnt > 40) {
			const Real occerr = ABS(occupation_err);
			if (occerr < 1.E-8) {
				WRN("The iteration stoped since the occerr is less than 1.E-8.");
				return true;
			}
		}
		return false;
}

bool NORG::green_modify_converged() const
{
	const Real cverr = ABS(correctionerr);
	if (cverr < 1.E-5) { return true; }
	return false;
}

bool NORG::green_modify_converged(Real correctionerr_i) const
{
	const Real cverr = ABS(correctionerr_i);
	if (cverr < 1.E-5) { return true; }
	return false;
}

//------------------------------------------------------------------ spin-spin ------------------------------------------------------------------

Real NORG::sz_imp_sz_bath(const Int imp_postition, const VecReal& vgs_i)
{
	VecPartition row_H(mm.np(), mm.id(), scsp.dim);
	Real sz_imp_sz_bath_i(0.);
	// VecReal gstate_part = vgs_i.truncate(row_H.bgn(), row_H.end());
	for_Int(h_i, row_H.bgn(), row_H.end())
	{
		StateStatistics a(h_i, scsp.wherein_NocSpace(h_i), scsp);
		Int conuter = 0;
		// for imp up spin:
#define ndiv scsp.ndivs
		if(a.occ_n[imp_postition][0] == 1) {
			conuter += SUM(a.occ_n[imp_postition].truncate(1,5));
			conuter -= SUM(a.occ_n[imp_postition + 1].truncate(1,5));

		}
		// for imp down spin:
		if(a.occ_n[imp_postition + 1][0] == 1){
			conuter -= SUM(a.occ_n[imp_postition].truncate(1,5));
			conuter += SUM(a.occ_n[imp_postition + 1].truncate(1,5));
		}
		sz_imp_sz_bath_i += conuter * vgs_i[h_i] * vgs_i[h_i];
	}
	Real sz_imp_sz_bath = mm.Allreduce(sz_imp_sz_bath_i) / 4.;
	return sz_imp_sz_bath;
}

//------------------------------------------------------------------ print out ------------------------------------------------------------------

void NORG::write_norg_info(Int iter_cnt) const {
	using namespace std;
	Real sc = 0;
	for_Int(i, 0, final_ground_state.size()) sc -= SQR(final_ground_state[i]) * log(SQR(final_ground_state[i]));
	OFS ofs_app_scorrelation(tox + "slater_correlation.txt", std::ios::app);
	ofs_app_scorrelation << iofmt("sci");
	ofs_app_scorrelation << setw(4) << iter_cnt;
	ofs_app_scorrelation << "\t" << setw(w_Real) << sc;
	ofs_app_scorrelation << endl; ofs_app_scorrelation.close();

	OFS ofs_app_energy(tox + "ground_energy.txt", std::ios::app);
	ofs_app_energy << iofmt("sci");
	ofs_app_energy << setw(4) << iter_cnt;
	ofs_app_energy << "\t" << setw(w_Real) << groune_lst;
	ofs_app_energy << endl; ofs_app_energy.close();

	// Only show for the upper spin.
	if(iter_cnt != 0)for_Int(i, 0, p.nband) {
		OFS ofs_app_occupation(tox +"imp"+ STR(i*2) + "_bath_occupation_number.txt", std::ios::app);
		ofs_app_occupation << iofmt("sci");
		ofs_app_occupation << setw(4) << iter_cnt;
		for_Int(j, 0, p.nI2B[i*2]) { ofs_app_occupation << "  " << setw(w_Real) << occnum_lst[SUM_0toX(p.nI2B, i*2) + j]; }
		ofs_app_occupation << endl;
		ofs_app_occupation.close();
	}
	
	// Only show for the upper spin.
	if(iter_cnt != 0)for_Int(i, 0, p.nband) {
		VEC<VecReal> dmeigen = oneedm.check_dm_get_occupation_number();
		OFS ofs_app_occupation_all(tox +"imp"+ STR(i*2) + "_occupation_number.txt", std::ios::app);
		ofs_app_occupation_all << iofmt("sci");
		ofs_app_occupation_all << setw(4) << iter_cnt;
		for_Int(j, 0, dmeigen[i].size()) { ofs_app_occupation_all << "  " << setw(w_Real) << dmeigen[i][j]; }
		ofs_app_occupation_all << endl;
		ofs_app_occupation_all.close();
	}

}


void NORG::write_state_info(Int iter_cnt) const {
	using namespace std;
	OFS ofs_app_state(tox + STR(iter_cnt) + "write_state_info.txt", ios::app);
	ofs_app_state << setw(4) << "iter_cnt";
	ofs_app_state << "\t" << setw(w_Real) << "state";
	ofs_app_state << "\t" << setw(w_Real) << "state_norm";
	ofs_app_state << "\t" << setw(w_Real) << "norm__sort";
	ofs_app_state << endl; 
	VecReal temp_norm_state = SQR(final_ground_state);
	VecReal temp_norm__sort = temp_norm_state; sort(temp_norm__sort);
	for_Int(i, 0, final_ground_state.size()) 
	{
		ofs_app_state << iofmt("sci");
		ofs_app_state << setw(4) << iter_cnt;
		ofs_app_state << "\t" << setw(w_Real) << final_ground_state[i];
		ofs_app_state << "\t" << setw(w_Real) << temp_norm_state[i];
		ofs_app_state << "\t" << setw(w_Real) << temp_norm__sort[final_ground_state.size()-1-i];
		ofs_app_state << endl; 
	}
	ofs_app_state.close();

}


Real NORG::return_nointeractions_ground_state_energy(const MatReal& h0_i) const {
	MatReal htemp = h0_i;
	VecReal vtemp(htemp.ncols(),0.);
	Real energy = 0.;
	heevr(htemp, vtemp);
	for_Int(i, 0, vtemp.size()) {
		if(vtemp[i]<0) energy += vtemp[i];
	}
	return energy;
}

//------------------------------------------------------------------- private -------------------------------------------------------------------

void NORG::show_the_nozero_number_of_tabel() 
{
	// Real size_one(oneedm.table[2].size()), size_two(n1mone.table[2].size()), size_tree(n1pone.table[2].size());
	// LLInt size_of_main_t(mm.Allreduce(size_one)), size_of_main_n1mt(mm.Allreduce(size_two)), size_of_main_n1pt(mm.Allreduce(size_tree));
	// if(mm) PIO(NAV3(size_of_main_t, size_of_main_n1mt, size_of_main_n1pt));
	Real size_one(oneedm.table[2].size());
	LLInt size_of_main_t(mm.Allreduce(size_one));
	if(mm) PIO(NAV(size_of_main_t));
}

MatReal NORG::save_transform_uormat(){
	MatReal transform_uormat(dmat(p.norbit, 1.));
	Int counter(0);
	if (!uormat.empty())
	{
		for (const auto& uormat_ii : uormat)
		{
			counter++;
			for_Int(i, 0, uormat_ii.nrows()) {
				for_Int(j, 0, uormat_ii.ncols()) {
					transform_uormat[i + counter][j + counter] = uormat_ii[i][j];
				}
			}
			counter += uormat_ii.nrows();
		}
	}
	return transform_uormat;
}
