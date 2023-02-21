/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

#include "operator.h"

using namespace std;


Operator::Operator(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& s_i):
	mm(mm_i), p(prmtr_i), scsp(s_i), table(find_h_idx())
{
}

Operator::Operator(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& s_i,const Tab &per_table):
	mm(mm_i), p(prmtr_i), scsp(s_i), table(per_table)
{
}

Operator::Operator(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& s_i, Str tab_name):
	mm(mm_i), p(prmtr_i), scsp(s_i), table(read_the_Tab(tab_name))
{
}

// (Deactivate) SparseH;
SparseMatReal Operator::find_hmlt()
{
	clock_t t_find_hmlt;
	TIME_BGN("find_hmlt" + NAV(mm.id()), t_find_hmlt);
	VecPartition row_H(mm.np(), mm.id(), scsp.dim);
	SparseMatReal hmlt_splited(row_H.len(), scsp.dim, mm);
	// if (mm) WRN("H begin and end" + NAV3(row_H.len(), row_H.bgn(), row_H.end()));
	
	Int numbercount(0);		//test variable
	for_Int(h_i, row_H.bgn(), row_H.end())
	{
		map<Idx, Real> va;
		StateStatistics a(h_i, scsp.wherein_NocSpace(h_i), scsp);
		// for the H_0 and H_imp (two fermi term diagonal) 
		va[h_i] = 0.;
		for (const auto &x : a.filled_spinless)
		{
			for(const auto &i : x)
			{
				va[h_i] += scsp.hopint[i][i];
			}
		}

#define Nup a.cfg.cf
#define Ndw a.cfg.cf
#define ndiv scsp.ndivs
		for_Int(i, 0, p.norg_sets)
		{
			va[h_i] -= scsp.mu[i] * Nup[i * ndiv][i] + scsp.mu[i] * Ndw[i * ndiv + 1][i];
			va[h_i] += scsp.u_hbd * Nup[i * ndiv][i] * Ndw[i * ndiv + 1][i];
			for_Int(j, 0, p.norg_sets) if (i != j)
			{
				va[h_i] += (scsp.u_hbd - 2 * scsp.j_ob) * Nup[i * ndiv][i] * Ndw[i * ndiv + 1][j];
				va[h_i] += (scsp.u_hbd - 3 * scsp.j_ob) * 0.5 * Nup[i * ndiv][i] * Nup[i * ndiv][j];
				va[h_i] += (scsp.u_hbd - 3 * scsp.j_ob) * 0.5 * Ndw[i * ndiv + 1][i] * Ndw[i * ndiv + 1][j];
			}
		}
		// off_diagonal_term
		// [i][0]:annihilation orbit's position; [i][1]:creation orbit's positon; [i][2]:Colum idx(i); [i][3]:sign(anticommutativity)
		vector<VecInt> off_dt;
		for_Int(idx_sets, 0, p.norg_sets){
			vector<VecInt> off_dt_next(a.find_each_spiless_group_off_diagonal_term(a.divocchop_ingroup(a.div_idx, idx_sets),idx_sets));
			off_dt.insert(off_dt.end(), off_dt_next.begin(), off_dt_next.end());
		} 
		for_Int(i, 0, off_dt.size()) {
			if (va.find(off_dt[i][2]) == va.end()) va[off_dt[i][2]] = 0.;
			else ERR("TO Think why off dianago has term?");
			va[off_dt[i][2]] += off_dt[i][3] * scsp.hopint[off_dt[i][1]][off_dt[i][0]];	// for the hopintmatrix C^+C at imp is 0.
		}
		Int count(0);
		for (auto it : va) if (it.second != 0.) { 
			hmlt_splited.addelement(it.second, it.first, h_i - row_H.bgn()); 
			++count;
		}
		if (count == 0)hmlt_splited.add_empty_row(h_i - row_H.bgn());
		numbercount += va.size();
	}
	TIME_END("find_hmlt" + NAV2(mm.id(),numbercount), t_find_hmlt);
	return hmlt_splited;
}

Tab Operator::find_h_idx()
{
	clock_t t_find_hmlt_table;
	t_find_hmlt_table = clock(); // TIME_BGN("find_hmlt_table" + NAV(mm.id()), t_find_hmlt_table);
	VecPartition row_H(mm.np(), mm.id(), scsp.dim);
	Tab h_idxs(3);
	MatInt mat_hop_pos(scsp.hopint.nrows(),scsp.hopint.ncols());
	for_Int(i, 0, mat_hop_pos.nrows()) for_Int(j, 0, mat_hop_pos.ncols()) mat_hop_pos[i][j] = i * mat_hop_pos.ncols() + j;
	Int h_hbd_idx(mat_hop_pos.size()	+ 1);
	Int h_orb_ud_idx(mat_hop_pos.size() + 2);
	Int h_orb_uu_idx(mat_hop_pos.size() + 3);
	Int h_orb_dd_idx(mat_hop_pos.size() + 4);
	for_Int(h_i, row_H.bgn(), row_H.end())
	{
		// To save as sparse matrix, [0]: row number;[1]: colum number;[2]: idx.
		VecInt h_idx(3, 0);
		Int sparse_idx(h_i - row_H.bgn());
		//WRN("wherein_NocSpace" + NAV(h_i - scsp.idx_div[scsp.wherein_NocSpace(h_i)]));
		StateStatistics a(h_i, scsp.wherein_NocSpace(h_i), scsp);
		for (const auto &x : a.filled_spinless) {
			for (const auto &i : x) {
				h_idx = { sparse_idx, h_i, mat_hop_pos[i][i] + 1 };
				for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
			}
		}


		#define ndiv scsp.ndivs
		{ // normal form of interation.
		// add the U.
		for_Int(i, 0, p.nband) if( a.cfg.cf[(i * 2) * ndiv].isocc(0) && a.cfg.cf[(i * 2 + 1) * ndiv].isocc(0) ){
			h_idx = { sparse_idx, h_i, h_hbd_idx};
			for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
		}

		// add the up-down term.
		for_Int(i, 0, p.nband) {
			for_Int(j, 0, p.nband) if(i != j && a.cfg.cf[(i * 2) * ndiv].isocc(0) && a.cfg.cf[(j * 2 + 1) * ndiv].isocc(0) ){
				h_idx = { sparse_idx, h_i, h_orb_ud_idx};
				for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
			}
		}

		// add the up-up term.
		for_Int(i, 0, p.nband) {
			for_Int(j, 0, p.nband) if(i != j && a.cfg.cf[(i * 2) * ndiv].isocc(0) && a.cfg.cf[(j * 2) * ndiv].isocc(0) ){
				h_idx = { sparse_idx, h_i, h_orb_uu_idx};
				for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
			}
		}

		// add the down-down term.
		for_Int(i, 0, p.nband) {
			for_Int(j, 0, p.nband) if(i != j && a.cfg.cf[(i * 2 + 1) * ndiv].isocc(0) && a.cfg.cf[(j * 2 + 1) * ndiv].isocc(0) ){
				h_idx = { sparse_idx, h_i, h_orb_dd_idx};
				for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
			}
		}
		}
		// // add the U.
		// for_Int(i, 0, p.nband) {
		// 	if((Real(a.cfg.cf[(i * 2) * ndiv][0]) - 0.5) * (Real(a.cfg.cf[(i * 2 + 1) * ndiv][0]) - 0.5) > 0 ) h_idx = { sparse_idx, h_i, h_hbd_idx };
		// 	else h_idx = { sparse_idx, h_i, -h_hbd_idx };
		// 	for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
		// }

		// // add the up-down term.
		// for_Int(i, 0, p.nband) {
		// 	for_Int(j, 0, p.nband) {
		// 		if (i != j && (Real(a.cfg.cf[(i * 2) * ndiv][0]) - 0.5) * (Real(a.cfg.cf[(j * 2 + 1) * ndiv][0]) - 0.5) > 0) h_idx = { sparse_idx, h_i, h_hbd_idx };
		// 		// if (i != j && a.cfg.cf[(i * 2) * ndiv].isocc(0) && a.cfg.cf[(j * 2 + 1) * ndiv].isocc(0)) h_idx = { sparse_idx, h_i, h_hbd_idx };
		// 		else h_idx = {sparse_idx, h_i, -h_orb_ud_idx};
		// 		for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
		// 	}
		// }

		// // add the up-up term.
		// for_Int(i, 0, p.nband) {
		// 	for_Int(j, 0, p.nband) {
		// 	if(i != j && (Real(a.cfg.cf[(i * 2) * ndiv][0]) - 0.5) * (Real(a.cfg.cf[(j * 2) * ndiv][0]) - 0.5) > 0) h_idx = { sparse_idx, h_i, h_orb_uu_idx};
		// 	else h_idx = { sparse_idx, h_i, -h_orb_uu_idx};
		// 		for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
		// 	}
		// }

		// // add the down-down term.
		// for_Int(i, 0, p.nband) {
		// 	for_Int(j, 0, p.nband) {
		// 	if(i != j && (Real(a.cfg.cf[(i * 2 + 1) * ndiv][0]) - 0.5) * (Real(a.cfg.cf[(j * 2 + 1) * ndiv][0]) - 0.5) > 0) h_idx = { sparse_idx, h_i, h_orb_dd_idx};
		// 	else h_idx = { sparse_idx, h_i, -h_orb_dd_idx};
		// 		for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
		// 	}
		// }
 
		// for_Int(i, 0, p.nband) {
		// 	if((Real(a.cfg.cf[0][i]) - 0.5) * (Real(a.cfg.cf[0][p.nband + i]) - 0.5) > 0 ) h_idx = { sparse_idx, h_i, h_hbd_idx };
		// 	else h_idx = { sparse_idx, h_i, -h_hbd_idx };
		// 	for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
		// }

/* 
		#define ndiv scsp.ndivs

		for_Int(i, 0, p.norg_sets)
		{
			for_Int(j, i, p.norg_sets)
			{
				if (i != j)
				{
					if (i / 2 == j / 2 && i + 1 == j)
					{
						if (a.cfg.cf[i * ndiv][0] && a.cfg.cf[j * ndiv][0])
							h_idx = {sparse_idx, h_i, h_hbd_idx};
						// if (a.cfg.cf[i * ndiv][0] != a.cfg.cf[j * ndiv][0])
						// 	h_idx = {sparse_idx, h_i, -h_hbd_idx};
						for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
					}
					if (i / 2 != j / 2 && i % 2 == 1 && j % 2 == 0)
					{
						if (a.cfg.cf[i * ndiv][0] && a.cfg.cf[j * ndiv][0])
							h_idx = {sparse_idx, h_i, h_orb_ud_idx};
						// if (a.cfg.cf[i * ndiv][0] != a.cfg.cf[j * ndiv][0])
						// 	h_idx = {sparse_idx, h_i, -h_orb_ud_idx};
						for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
					}
					if (i / 2 != j / 2 && i % 2 == 0 && j % 2 == 1)
					{
						if (a.cfg.cf[i * ndiv][0] && a.cfg.cf[j * ndiv][0])
							h_idx = {sparse_idx, h_i, h_orb_ud_idx};
						// if (a.cfg.cf[i * ndiv][0] != a.cfg.cf[j * ndiv][0])
						// 	h_idx = {sparse_idx, h_i, -h_orb_ud_idx};
						for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
					}
					if (i / 2 != j / 2 && i % 2 == 1 && j % 2 == 1)
					{
						if (a.cfg.cf[i * ndiv][0] && a.cfg.cf[j * ndiv][0])
							h_idx = {sparse_idx, h_i, h_orb_uu_idx};
						// if (a.cfg.cf[i * ndiv][0] != a.cfg.cf[j * ndiv][0])
						// 	h_idx = {sparse_idx, h_i, -h_orb_uu_idx};
						for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
					}
					if (i / 2 != j / 2 && i % 2 == 0 && j % 2 == 0)
					{
						if (a.cfg.cf[i * ndiv][0] && a.cfg.cf[j * ndiv][0])
							h_idx = {sparse_idx, h_i, h_orb_dd_idx};
						// if (a.cfg.cf[i * ndiv][0] != a.cfg.cf[j * ndiv][0])
						// 	h_idx = {sparse_idx, h_i, -h_orb_dd_idx};
						for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
					}
				}
			}
		}
 */

		for_Int(idx_sets, 0, p.norg_sets)
		{
			// off_diagonal_term
			// i[0]:annihilation orbit's position; i[1]:creation orbit's positon; i[2]:Colum idx(i); i[3]:sign(anticommutativity)
			VEC<VecInt> off_dt_next(a.find_each_spiless_group_off_diagonal_term(a.divocchop_ingroup(a.div_idx, idx_sets), idx_sets));
			for (const auto &i : off_dt_next){
				h_idx = {sparse_idx, i[2], i[3] * (mat_hop_pos[i[1]][i[0]] + 1)};
				for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
			}
		}
	}
	// TIME_END("t_find_hmlt_table" + NAV(mm.id()), t_find_hmlt_table);
	for_Int(jj, 0, 3) h_idxs[jj].shrink_to_fit();
	return std::move(h_idxs);
}

SparseMatReal Operator::find_hmlt(const Tab h_idx) const
{
	// clock_t t_find_hmlt;
	// TIME_BGN("find_hmlt" + NAV(mm.id()), t_find_hmlt);
	VecPartition row_H(mm.np(), mm.id(), scsp.dim);
	SparseMatReal hmlt_splited(row_H.len(), scsp.dim, mm);
	Real diagonal(0.);	Bool flag(false);
	for_Idx(pos, 0, h_idx[2].size())
	{
		Int coefficient_comm = h_idx[2][pos] >= 0 ? 1 : -1;
		if (h_idx[1][pos] == h_idx[0][pos] + row_H.bgn()) {diagonal += coefficient_comm * scsp.coefficient[abs(h_idx[2][pos])];flag=true;}
		else {
			if (flag) {hmlt_splited.addelement(diagonal, h_idx[1][pos-1], h_idx[0][pos-1]); diagonal = 0.;flag=false;}
			hmlt_splited.addelement(coefficient_comm * scsp.coefficient[abs(h_idx[2][pos])], h_idx[1][pos], h_idx[0][pos]);
		}
	}
	// TIME_END("find_hmlt" + NAV(mm.id()), t_find_hmlt);
	hmlt_splited.shrink_to_fit();
	return hmlt_splited;
}


// Try to find the Ground state.
MatReal Operator::lowest_eigpairs(const Idx n, bool if_need_fast, Int wish_nev)
//	n      ( input) is the dimension of the matrix
//	wish_nev    ( input) is the expected maximum eigenstates to be found
//	actual number of eigenpairs is returned
//	nev_eff is the number of eigenvectors corresponding to complete degenerate levels reached
{
	UURandomSFMT uur;
	VecReal eval(wish_nev, 0.);				// is the eigenvalues
	// VEC < VecReal> eigenvec_i;			// a VEC to contain all the eigenvector.
	MatReal eigenvec_i(wish_nev, n);		// a VEC to contain all the eigenvector.
	VecInt ev_dgcy(wish_nev, 0.);			// cotain all the degeneracy eigenvector's number for each eigenvalue.
	// if (mm)PIO("find_hmlt BEGIN :::");
	SparseMatReal sep_hmltoperator = find_hmlt(table);
	// if(mm) WRN("finished fiding the hmlt")
//#ifdef _ASSERTION_
//	if (mm.np() == 1) {
//		WRN("test begin:::");
//		if (sep_hmltoperator.if_hermitian() == false)ERR("This sparse matrix is not Hermitian matrix!!!!");
//		WRN("test end  !!!");
//	}
//#endif
	VecReal inital_state(n, 0.); uur(inital_state); inital_state -= 0.5;
	VecInt krylov_space_size = lanczos(eval, eigenvec_i, ev_dgcy, n, wish_nev, sep_hmltoperator, inital_state, if_need_fast, 9999);
	if(mm) {cout <<"PIO: krylov_space_size = "; for_Int(i, 0, krylov_space_size.size()) cout << krylov_space_size[i] << "; "; cout << std::endl;}
	// if(mm) std::cout << "The eigenvalue" << iofmt("sci") << eval << std::endl;
	groundstate_energy = eval[0];
	// MatReal eigenvec(eigenvec_i.size(),n);
	// for_Int(i, 0, eigenvec_i.size()) eigenvec[i] = eigenvec_i[i];
	ground_state = eigenvec_i[0];
	// if(mm) WRN("finished lanczos.")
#ifdef _ASSERTION_
	//WRN("TEST_Lanczos is right:::" + NAV(mm.np()));
	//VecReal test_a(eigenvec[0]);
	//VecReal test_b(sep_hmltoperator * test_a);
	//WRN("TEST_Lanczos[0] state" + NAV3(test_a.isnormalized(), eval[0], test_b.avg_abs_elem_diff(eval[0] * test_a)));
#endif
	// if(mm) WRN(NAV(groundstate_energy));
	return eigenvec_i;
}


VecReal Operator::sn_prtcl_ex_state(const Int imp_div, const VecReal ground_state, const Int crtann) const
{
	VecPartition row_H(mm.np(), mm.id(), scsp.dim);
	VecReal ex_state_part(row_H.len(), 0.);
	for_Int(h_i, row_H.bgn(), row_H.end())
	{
		Int subscsp(scsp.wherein_NocSpace(h_i));
		ComDivs cfg(h_i - scsp.idx_div[subscsp], (scsp.div[subscsp]), (scsp.sit_mat), true);
		if (crtann == +1)if (cfg.cf[imp_div * scsp.ndivs].isuno(0))ex_state_part[h_i - row_H.bgn()] = 1.;
		if (crtann == -1)if (cfg.cf[imp_div * scsp.ndivs].isocc(0))ex_state_part[h_i - row_H.bgn()] = 1.;
	}
	ex_state_part *= ground_state.truncate(row_H.bgn(), row_H.end());
	VecReal ex_state(scsp.dim, 0.);
	ex_state = mm.Allgatherv(ex_state_part, row_H);
	return std::move(ex_state);
}

// (deprecated!!!)ONLY use for the imp_postition.
VecReal Operator::particle_number_Inner_product(const Int imp_div, const Int crtann) const
{
	VecPartition row_H(mm.np(), mm.id(), scsp.dim);
	VecReal ex_state_part(row_H.len(), 0.);
	for_Int(h_i, row_H.bgn(), row_H.end())
	{
		Int subscsp(scsp.wherein_NocSpace(h_i));
		ComDivs cfg(h_i - scsp.idx_div[subscsp], (scsp.div[subscsp]), (scsp.sit_mat), true);
		if (crtann == +1)if (cfg.cf[imp_div].isuno(0))ex_state_part[h_i - row_H.bgn()] = 1.;
		if (crtann == -1)if (cfg.cf[imp_div].isocc(0))ex_state_part[h_i - row_H.bgn()] = 1.;
	}
	VecReal ex_state(scsp.dim, 0.);
	ex_state = mm.Allgatherv(ex_state_part, row_H);
	return ex_state;
}

//------------------------------------------------------------------ io ------------------------------------------------------------------

void Operator::save_the_Tab(Tab& tab, Str name) const{
	Int size_temp(tab[0].size());
	Int size = mm.Allreduce(size_temp);
	VecInt v_size_i(1); v_size_i[0] = tab[0].size();
	VecPartition split_v_size(mm.np(), mm.id(), mm.np());
	VecInt v_size = mm.Allgatherv(v_size_i, split_v_size);
	if(mm) {// write the Tab's size's info
		OFS ofs;	ofs.open(name + ".inf");
		ofs << setw(9) << "dim" << setw(p_Real) << scsp.dim << endl;
		ofs << setw(9) << "size" << setw(p_Real) << size << endl;
		for_Int(i, 0, mm.np())	{
			ofs << setw(9) << "size_np"+STR(i) << setw(p_Real) << v_size[i] << endl;
		}
		ofs.close();
	}

	OFS ofs;	ofs.open(name + ".bdat");
	for_Int(i, 0, tab.size()) {
		VecPartition split_table(mm.np(), mm.id(), size, v_size);
		VecInt temp_tabi = mm.Gatherv(Vec(tab[i]), split_table);
		if(mm) {
			biwrite(ofs, CharP(temp_tabi.p()), temp_tabi.szof());
			
			// if(mm) WRN(NAV2(temp_tabi.size(),temp_tabi.truncate(size-100,size)));
		}		
	}
	ofs.close();
}

Tab Operator::read_the_Tab(Str name) const{
	
	// WRN("Here is fine"+NAV(name));
	VecInt v_size(mm.np(), 0);
	Int size(-1);	Tab tab(3);
	{
		IFS ifs(STR(name + ".inf"));	Str strr;
		while(1) {// read the Tab's size's info
			ifs >> strr;
			if(strr == "size")	ifs >> size;
			for_Int(i, 0, mm.np()) if(strr == "size_np"+STR(i))	ifs >> v_size[i];
			if (!ifs) break;
		}
	}
	// if(mm) WRN("Here is fine"+NAV(size));
	{
		IFS ifs(STR(name + ".bdat"));
		// VecPartition split_tab(mm.np(), mm.id(), size);		
		VecPartition split_tab(mm.np(), mm.id(), size, v_size);
		for_Int(i, 0, tab.size()) {
			VecInt temp(size); biread(ifs, CharP(temp.p()), temp.szof());
			// if(mm) WRN(NAV(temp.truncate(size-100,size)));
			tab[i] = temp.truncate(split_tab.bgn(),split_tab.end()).stdvec();
			// SWAP(tab[i], temp.truncate(split_tab.bgn(),split_tab.end()).stdvec());
		}
	}
	return tab;
}