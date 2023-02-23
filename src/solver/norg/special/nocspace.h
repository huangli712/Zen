#pragma once

/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

#include "specs.h"
#include "prmtr.h"
#include "impurity.h"
// impurity model
// At the Shortcut space(NocSpace) we set first divison as impurity, and the active orbital(bath site) on the middle.  

class NocSpace {

private:
	const BnmlFast bf;
public:
	const MyMpi &mm;				// parameters
	const Prmtr& p;					// parameters
	Idx ndivs;						// The amount of divisons's number. 
	MatInt control_divs;			// to set the number of division and the shortcut restratin.
	VecInt shortcut_countvec;		// the SUM in row for each colum.
	MatInt sit_mat;					// spinless - orbits number in each division.
	Idx nspa;						// The amount of partical's number.
	VecInt nppso;					// nppso mean: number of partical per spin orbital.
	
	Idx dim;						// the dimension of the Shortcut space.
	VEC<MatInt> div;				// a set of combined number.
	std::map<std::string, Idx> divs_to_idx;	// a set of combined number.
	VEC<Int> idx_div;				// a set of idx(The idx of the begining with 0) for each subspace.
	VecReal mu;						// -\mu_{k}^{\prime}b_{k,\sigma}^{+}b_{k,\sigma}
	VEC<VecReal> t_ose;				// H_0 onset energy
	VEC<VecReal> t_hyb;				// hopping parameter from bath to imp.
	MatReal hopint;					// transformed hopping integral
	VecReal coefficient;			// coefficient for all the H's terms, consturct by t_ose and t_hyb.
	
	
	Real u_hbd;						// The Hubbard term U.
	Real j_ob;						// The Hunter couping term J.
	
private:
	void set_control();
	// Find all the combined number subspaces.
	void find_combined_number_subspaces(const Int mode = 0);
	// // Find all the combined number subspaces by nppso.
	// void find_combined_number_subspaces(const VecInt& nppso);

	// Find all the combined number subspaces, with speed up.
	void find_all_noc_subspaces();
	void find_all_noc_subspaces_by_row();

	void print(std::ostream& os, const Str& var, const Str& val, const Str& comment) const {
		using namespace std;
		Str true_var = var == "\"\"" ? "" : var;
		os << rght_justify(true_var, 16) << (true_var == "" ? "   " : " = ") << left_justify(val, w_Real)
			<< "    # " + comment << endl;
	}
	bool if_div_in_restraint(const VecInt& restraint, const Int position, const Int max, const Int now) const;
	
	bool if_col_divs_in_restraint(const Int& restraint, const VEC<Int>& divcol_i, Int col_idx) const;
	bool if_row_divs_in_restraint(const Int& restraint, const VEC<Int>& divrow_i, VecInt count_sit_mat) const;

	bool ifin_NocSpace_judge_by_nppso(const MatInt& spilss_div, const VecInt& nppso) const;

	VecInt multi_judger(const VEC<VEC<int> >& s, const VEC<VEC<int> >& a) const;
	VecInt multi_judger_by_row(const VEC<VEC<int> >& s, const VEC<VEC<int> >& a) const;

	bool check_each_column(const Int& col_pos, const VecInt& div_colsum) const {
		if (col_pos == div_colsum.size() / 2) return true;
		if (col_pos < div_colsum.size() / 2 && div_colsum[col_pos] >= shortcut_countvec[col_pos] + control_divs[0][col_pos] && div_colsum[col_pos] <= shortcut_countvec[col_pos]) return true;
		if (col_pos > div_colsum.size() / 2 && div_colsum[col_pos] >= 0 && div_colsum[col_pos] <= control_divs[0][col_pos]) return true;
		// if (control_divs[0][col_pos] ==0 && ((col_pos < div_colsum.size()/2 && div_colsum[col_pos]==shortcut_countvec[col_pos]) || (col_pos > div_colsum.size()/2 && div_colsum[col_pos]==0))) return true;
		return false;
	}

	bool check_correlated_column(const Int& col_pos, const VecInt& div_colsum) const {
		if (shortcut_countvec[col_pos] - div_colsum[col_pos] + div_colsum[div_colsum.size() - col_pos] <= control_divs[0][div_colsum.size() - col_pos]) return false;
		else return true;
	}

	VecInt read_from_col_lable(const VEC<Int> x, const VEC<VEC<Int> > a) const;

	void find_all_possible_state(VEC<VEC<Int> >& a, VEC<VEC<Int> >& s) const;
	void find_all_possible_state_by_row(VEC<VEC<Int> >& a, VEC<VEC<Int> >& s) const;

	VEC<VEC<Int> > cart_product(const VEC<VEC<int> >& v)const;

	VEC<VEC<Int> > cart_product_monitor_col(const VEC<VEC<int> >& v, const VEC<VEC<Int> >& a)const;
	VEC<VEC<Int> > cart_product_monitor_row(const VEC<VEC<int> >& v, const VEC<VEC<int> >& a) const;

	Idx read_the_Tab(Str name) const{
		Idx temp_dim(-1);	
		IFS ifs(STR(name + ".inf"));	Str strr;
		while(1) {// read the Tab's size's info
			ifs >> strr;
			if(strr == "dim")	ifs >> temp_dim;
			if (!ifs) break;
		}
		// WRN(NAV(temp_dim))
		return temp_dim;
	}
	
	Str nppso_str(const VecInt &nppso_i) const	{
		Str temp; for_Int(i, 0, nppso_i.size()) if (i % 2 == 0) temp += "-" + STR(nppso_i[i]);	return temp;
	}
public:
	// It assume that we already have the hopint from the Impurity class, but still have not rotated it yet.
	VecReal set_row_primeter_by_gived_mat(const VEC<MatReal>& uormat_i, const MatReal& h0);
	// Expand the Shortcut space under Number of particles(NumberSpa).
	NocSpace(const MyMpi& mm_i, const Prmtr& prmtr_i, const Int& NumberSpa);
	NocSpace(const MyMpi& mm_i, const Prmtr& prmtr_i, const MatReal& imp_i_h0, const Int& NumberSpa);
	NocSpace(const MyMpi& mm_i, const Prmtr& prmtr_i, const MatReal& imp_i_h0, const VecInt& nppso_i);
	NocSpace(const MyMpi& mm_i, const Prmtr& prmtr_i, const MatReal& imp_i_h0, const VecInt& nppso_i, Str tab_name);
	
	// bool ifin_NocSpace(VecInt& ud) const;
	bool ifin_NocSpace(MatInt& ud) const;
	// bool ifin_NocSpace_for_green(MatInt& spilss_div, const VecInt& nppso, const Int& crtorann) const;
	bool ifin_NocSpace(MatInt& spilss_div, const VecInt& nppso) const;
	
	bool ifin_NocSpace_more_strict(MatInt& spilss_div, const VecInt& nppso) const;

	Int wherein_NocSpace(const Int& h_i) const;

	void print(std::ostream& os = std::cout) const;


	VecInt free_div_base_decode(Int idx, VEC<VEC<Int> > v) const;
};
