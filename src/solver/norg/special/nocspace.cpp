/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

#include "nocspace.h"
#define nip p.norg_sets
#define i2b p.nI2B
#define nob p.norbit
// the NocSpace, is start frome the shortcutspace.
NocSpace::NocSpace(const Prmtr& prmtr_i, const Int& NumberSpa) :
	p(prmtr_i), ndivs(prmtr_i.ndiv), nspa(NumberSpa),sit_mat(prmtr_i.control_divs.truncate_row(1, control_divs.nrows())),
	hopint(nob, nob, 0.), coefficient((p.norbit * p.norbit) + 5, 0.),
	control_divs(p.norg_sets + 1, prmtr_i.ndiv, 0), shortcut_countvec(prmtr_i.ndiv, 0), dim(0)
{
	set_control();
	find_combined_number_subspaces();
	//WRN("Finish finding the subspace and the Dim is:" + NAV(dim));
}

NocSpace::NocSpace(const Prmtr& prmtr_i, const MatReal& imp_i_h0, const Int& NumberSpa) :
 	p(prmtr_i), ndivs(prmtr_i.ndiv), nspa(NumberSpa),sit_mat(prmtr_i.control_divs.truncate_row(1, control_divs.nrows())),
 	hopint(imp_i_h0), coefficient((p.norbit * p.norbit) + 5, 0.),
 	control_divs(p.norg_sets + 1, prmtr_i.ndiv, 0), shortcut_countvec(prmtr_i.ndiv, 0), dim(0)
{
 	set_control();
 	find_combined_number_subspaces();
}

// nppso mean: number of partical per spin orbital.
NocSpace::NocSpace(const Prmtr& prmtr_i, const MatReal& imp_i_h0, const VecInt& nppso_i) :
	p(prmtr_i), ndivs(prmtr_i.ndiv), nspa(SUM(nppso_i)),sit_mat(prmtr_i.control_divs.truncate_row(1, control_divs.nrows())),
	hopint(imp_i_h0), coefficient((p.norbit * p.norbit) + 5, 0.), nppso(nppso_i),
	control_divs(p.norg_sets + 1, prmtr_i.ndiv, 0), shortcut_countvec(prmtr_i.ndiv, 0), dim(0)
{
	set_control();
	// find_combined_number_subspaces(1);
	find_all_noc_subspaces();
	WRN("Begin find_combined_number_subspaces()"+ NAV(dim));
}

void NocSpace::set_control()
{
	control_divs = p.control_divs;
	for_Int(j, 0, control_divs.ncols())for_Int(i, 1, control_divs.nrows()) shortcut_countvec[j] += control_divs[i][j];
	// sit_mat.reset(control_divs.truncate_row(1, control_divs.nrows()));
}


VecReal NocSpace::set_row_primeter_by_gived_mat(const VEC<MatReal>& uormat_i, const MatReal& h0)
{
	hopint = h0;
	MatReal transform_uormat(dmat(p.norbit, 1.));
	Int counter(0);
	if (!uormat_i.empty())
	{
		for (const auto& uormat_ii : uormat_i)
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
	hopint = transform_uormat * hopint * transform_uormat.ct();
	VecReal coefficient_i(hopint.vec());
	VecReal check_point(1, 0.); coefficient_i.reset(concat(check_point, coefficient_i));
	// // VecReal orbital_correction({ u_hbd, (u_hbd - (2 * j_ob)), (u_hbd - (3 * j_ob)), (u_hbd - (3 * j_ob)) });
	// // VecReal orbital_correction({ u_hbd, p.hubbU12, p.hubbU12, p.hubbU12 });
	// {U, U':u d, U':u u, U':d d}
	VecReal orbital_correction({ p.hubbU, p.uprime, p.uprime - p.jz, p.uprime - p.jz });
	// VecReal orbital_correction({ 0.25 * p.hubbU, 0. * p.hubbU, 0. * p.hubbU, 0. * p.hubbU });
	return concat(coefficient_i, orbital_correction);
}

// find all the combined number subspaces OR find all the combined number subspaces with special partical control.
void NocSpace::find_combined_number_subspaces(const Int mode)
{
	Idx length_ECNS;			// length for each combined number subspaces(ECNS).
	idx_div.push_back(dim);
	VEC<VEC<Int> > a;
	for_Int(row, 1, control_divs.nrows()) {
		for_Int(col, 0, control_divs.ncols()) {
			VEC<Int> one_div;
			for_Int(spl, 0, control_divs[row][col] + 1) {
				if(if_div_in_restraint(control_divs[0],col, control_divs[row][col],spl)) one_div.push_back(spl);
			}
			a.push_back(one_div);
		}
	}
	Idx total_posibile(1);
	for (const auto& u : a) total_posibile *= u.size();

	if(mode == 1) {
		for_Idx(x, 0, total_posibile) {
			VecInt spilss_div_v(free_div_base_decode(x, a));
			MatInt spilss_div(spilss_div_v.mat(control_divs.nrows() - 1, control_divs.ncols()));
			if (ifin_NocSpace(spilss_div, nppso)) {
				div.push_back(spilss_div);	length_ECNS = 1;
				for_Int(i, 1, control_divs.nrows()) for_Int(j, 0, ndivs) length_ECNS = length_ECNS * bf(control_divs[i][j], spilss_div[i - 1][j]);
				dim += length_ECNS;			idx_div.push_back(dim);
			}
		}
	}

	if(mode == 0) {
		for_Idx(x, 0, total_posibile) {
			VecInt spilss_div_v(free_div_base_decode(x, a));
			MatInt spilss_div(spilss_div_v.mat(control_divs.nrows() - 1, control_divs.ncols()));
			if (ifin_NocSpace(spilss_div)) {
				div.push_back(spilss_div);	length_ECNS = 1;
				for_Int(i, 1, control_divs.nrows()) for_Int(j, 0, ndivs) length_ECNS = length_ECNS * bf(control_divs[i][j], spilss_div[i - 1][j]);
				dim += length_ECNS;			idx_div.push_back(dim);
			}
		}
	}
}

void NocSpace::find_all_noc_subspaces()
{
	Idx length_ECNS;						// length for each combined number subspaces(ECNS).
	idx_div.push_back(dim);
	VEC<VEC<Int> > a;
	VEC<VEC<Int> > s;

	find_all_possible_state(a, s);
	for (const auto &x : s)
	{
		VecInt spilss_div_v(read_from_col_lable(x,a));
		MatInt spilss_div_tr(spilss_div_v.mat(control_divs.ncols(), control_divs.nrows() - 1));
		MatInt spilss_div(spilss_div_tr.tr());
		if (ifin_NocSpace(spilss_div, nppso)) {
			div.push_back(spilss_div);
			divs_to_idx.insert(std::pair<std::string, Int>(spilss_div.vec().string(), dim));
			length_ECNS = 1;
			for_Int(i, 1, control_divs.nrows()) for_Int(j, 0, ndivs) length_ECNS = length_ECNS * bf(control_divs[i][j], spilss_div[i - 1][j]);
			dim += length_ECNS;
			idx_div.push_back(dim);
		}
	}
}

VecInt NocSpace::read_from_col_lable(const VEC<Int> x, const VEC<VEC<Int> > a) const{
	VecInt ren;
	for (const auto &lable : x)
	{
		VecInt ren_i(ren);
		VecInt col(a[lable]);
		ren.reset(concat(ren_i, col));
	}
	return ren;
}

void NocSpace::find_all_possible_state(VEC<VEC<Int> >& a, VEC<VEC<Int> >& s) const
{
	VEC<VEC<Int> > a_lable;
	Int counter(0);
	for_Int(col, 0, control_divs.ncols())
	{
		VEC<VEC<Int>> temp_a;
		VEC<VEC<Int> > a_rol_temp;
		for_Int(row, 1, control_divs.nrows())
		{
			VEC<Int> one_div;
			for_Int(spl, 0, control_divs[row][col] + 1)
			{
				if (if_div_in_restraint(control_divs[0], col, control_divs[row][col], spl))
					one_div.push_back(spl);
			}
			a_rol_temp.push_back(one_div);
		}
		temp_a = cart_product(a_rol_temp);
		for (const auto &one_divs : temp_a)	{
			if (if_col_divs_in_restraint(control_divs[0][col], one_divs, col))
				a.push_back(one_divs);
		}
		VEC<Int> a_lable_i;
		for_Int(i, counter, a.size()){
			a_lable_i.push_back(i);
			counter++;
		}
		a_lable.push_back(a_lable_i);
		//  WRN(NAV2( col,a.size()));
	}

	s = cart_product(a_lable);
}

bool NocSpace::ifin_NocSpace(MatInt& spilss_div) const
{
	if (SUM(spilss_div) == nspa) {
		for_Int(i, 1, control_divs.nrows()) for_Int(j, 0, ndivs) if (spilss_div[i - 1][j] < 0 || spilss_div[i - 1][j] > control_divs[i][j]) return false;
		VecInt spilsdiv_countvec(ndivs, 0);
		for_Int(j, 0, spilss_div.ncols())for_Int(i, 0, spilss_div.nrows()) spilsdiv_countvec[j] += spilss_div[i][j];
		if (spilsdiv_countvec[1] >= shortcut_countvec[1] + control_divs[0][1] && \
			spilsdiv_countvec[1] <= shortcut_countvec[1] && \
			spilsdiv_countvec[2] >= shortcut_countvec[2] + control_divs[0][2] && \
			spilsdiv_countvec[2] <= shortcut_countvec[2] && \
			spilsdiv_countvec[4] >= 0 && \
			spilsdiv_countvec[4] <= control_divs[0][4] && \
			spilsdiv_countvec[5] >= 0 && \
			spilsdiv_countvec[5] <= control_divs[0][5]) return true;
	}
	return false;
}

bool NocSpace::ifin_NocSpace(MatInt& spilss_div, const VecInt& nppso) const
{
	if(ndivs%2) ERR("ndivs is not a even number!");
	for_Int(i, 0, p.norg_sets)	if(SUM(spilss_div[i]) != nppso[i]) return false;
	if (SUM(spilss_div) == nspa) {
		for_Int(i, 1, control_divs.nrows()) for_Int(j, 0, ndivs) if (spilss_div[i - 1][j] < 0 || spilss_div[i - 1][j] > control_divs[i][j]) return false;
		VecInt spilsdiv_countvec(ndivs, 0);
		for_Int(j, 0, spilss_div.ncols())for_Int(i, 0, spilss_div.nrows()) spilsdiv_countvec[j] += spilss_div[i][j];
		/*
		if (spilsdiv_countvec[1] >= shortcut_countvec[1] + control_divs[0][1] && \
			spilsdiv_countvec[1] <= shortcut_countvec[1] && \
			spilsdiv_countvec[2] >= shortcut_countvec[2] + control_divs[0][2] && \
			spilsdiv_countvec[2] <= shortcut_countvec[2] && \
			spilsdiv_countvec[4] >= 0 && \
			spilsdiv_countvec[4] <= control_divs[0][4] && \
			spilsdiv_countvec[5] >= 0 && \
			spilsdiv_countvec[5] <= control_divs[0][5]){
				
				if(p.nooc_mode == STR("nooc"))return true;
				else if(p.nooc_mode == STR("cnooc")){
					if( shortcut_countvec[2] - spilsdiv_countvec[2] + spilsdiv_countvec[4] <= control_divs[0][4] && \
						shortcut_countvec[1] - spilsdiv_countvec[1] + spilsdiv_countvec[5] <= control_divs[0][5]) return true;
				}
				else ERR("nooc_mode in put was wrong!");
		}
		*/
		for_Int(i, 1, ndivs) if(!check_each_column(i, spilsdiv_countvec)) return false;
		if(p.nooc_mode == STR("nooc")) return true;
		else if(p.nooc_mode == STR("cpnooc")) {for_Int(i, 1, ndivs/2 - 1) if(check_correlated_column(i, spilsdiv_countvec)) return false;}
		else if(p.nooc_mode == STR("cnooc")) {for_Int(i, 1, ndivs/2) if(check_correlated_column(i, spilsdiv_countvec)) return false;}
		else ERR("nooc_mode in put was wrong!");
		return true;
	}
	else return false;
}

bool NocSpace::ifin_NocSpace_more_strict(MatInt& spilss_div, const VecInt& nppso) const
{
	for_Int(i, 0, p.norg_sets)	if(SUM(spilss_div[i]) != nppso[i]) return false;
	if (SUM(spilss_div) == nspa) {
		for_Int(i, 1, control_divs.nrows()) for_Int(j, 0, ndivs) if (spilss_div[i - 1][j] < 0 || spilss_div[i - 1][j] > control_divs[i][j]) return false;
		VecInt spilsdiv_countvec(ndivs, 0);
		for_Int(j, 0, spilss_div.ncols())for_Int(i, 0, spilss_div.nrows()) spilsdiv_countvec[j] += spilss_div[i][j];
		/*
		if (spilsdiv_countvec[1] >= shortcut_countvec[1] + control_divs[0][1] && \
			spilsdiv_countvec[1] <= shortcut_countvec[1] && \
			spilsdiv_countvec[2] >= shortcut_countvec[2] + control_divs[0][2] && \
			spilsdiv_countvec[2] <= shortcut_countvec[2] && \
			spilsdiv_countvec[4] >= 0 && \
			spilsdiv_countvec[4] <= control_divs[0][4] && \
			spilsdiv_countvec[5] >= 0 && \
			spilsdiv_countvec[5] <= control_divs[0][5] && \
			shortcut_countvec[1] - spilsdiv_countvec[1] + spilsdiv_countvec[5] <= control_divs[0][5]) return true;
		*/
		for_Int(i, 1, ndivs) if(check_each_column(i, spilsdiv_countvec)) return false;
	//  && shortcut_countvec[1] - spilsdiv_countvec[1] + spilsdiv_countvec[5] <= control_divs[0][5]
		return true;
	}
	else return false;
}

bool NocSpace::ifin_NocSpace_judge_by_nppso(const MatInt& spilss_div, const VecInt& nppso) const
{
	for_Int(i, 0, p.norg_sets)	if(SUM(spilss_div[i]) != nppso[i]) return false;

	return true;
}

bool NocSpace::if_div_in_restraint(const VecInt& restraint, const Int position, const Int max, const Int now) const
{
	if (position == 0 || position == ndivs / 2) return true;
	if (restraint[position] <= 0 && now >= (max + restraint[position])) return true;
	if (restraint[position] >= 0 && now <= restraint[position]) return true;
	return false;
}

bool NocSpace::if_col_divs_in_restraint(const Int& restraint, const VEC<Int>& divcol_i, Int col_idx) const
{
	VecInt divcol(divcol_i);
	// for_Int(i, 0, divcol.size()) if(divcol[i] < 0) return 0;
	if(restraint == 0) return true;
	if(restraint < 0 && SUM(divcol) >= (shortcut_countvec[col_idx] + restraint) ) return true;
	if(restraint > 0 && SUM(divcol) <= restraint ) return true;
	return false;
}

Int NocSpace::wherein_NocSpace(const Int& h_i)const
{
	Int comdiv_idx(0);
	for_Int(j, 0, idx_div.size()) {
		if (h_i < idx_div[j]) {
			comdiv_idx = j - 1;
			break;
		}
		else if (h_i >= idx_div[idx_div.size() - 1]) comdiv_idx = idx_div.size() - 1;
	}
	return comdiv_idx;
}


void NocSpace::print(std::ostream& os) const
{
#define nocspace_print(var, comment) print(os, NAME(var), STR(var), comment)

	using namespace std;
	Str cnooc = p.nooc_mode;


	os << "// NORG setting" << endl;

	// nocspace_print(ndivs, "The amount of divisons's number. ");
	nocspace_print(ndivs, "The amount of divisons's number. ");
	nocspace_print(cnooc, "Correlation nature orbital occupation constraint.");
	nocspace_print(control_divs, "to set the number of division and the shortcut restratin.");
	nocspace_print(nspa, "The amount of partical's number.");
	nocspace_print(dim, "the dimension of the Shortcut space.");
	// nocspace_print(h0, "transformed hopping integral");
	// nocspace_print(mu, "-mu");
	nocspace_print(p.hubbU, " The Hubbard term U.");
	// u_hbd, p.hubbU12, p.hubbU12, p.hubbU12

	os << "// prmtr print end  " << present() << endl;

#undef nocspace_print
}

VEC<VEC<Int> > NocSpace::cart_product (const VEC<VEC<Int> >& v)const
{
    VEC<VEC<Int> > s = {{}};
    for (const auto& u : v) {
        VEC<VEC<Int> > r;
        for (const auto& x : s) {
            for (const auto y : u) {
                r.push_back(x);
                r.back().push_back(y);
            }
        }
        s = move(r);
    }
    return s;
}

VecInt NocSpace::free_div_base_decode(Int idx, VEC<VEC<Int> > v) const
{
	VecInt base(v.size() + 1, 1);
	for_Int(i, 1, base.size()) base[i] = base[i - 1] * v[i - 1].size();

	VecInt rep(base.size() - 1, 0);
	for_Int(i, 0, rep.size()) rep[i] = v[i][(idx % base[i + 1]) / base[i]];
	if (idx >= base[base.size() - 1]) ERR(STR("free div base decode error."));
	return rep;
}