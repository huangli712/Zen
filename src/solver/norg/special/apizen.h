#pragma once

/*

coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2023
*/

#include "specs.h"
#include "prmtr.h"
#include "model.h"
#include "green.h"
#include "hyberr.h"
#include "bath.h"
#include "impurity.h"
#include "norg.h"

// The API for the ZEN framework as a pure impurity solver.

class APIzen {
public:
	const MyMpi& mm;
	const Prmtr& p;
	const Model& mdl;
	Int iter_cnt;
	VecReal se_err;
	
	Real p_n;				// partical number
	Real m_m;				// magnetic moment
	Real d_p;				// d-wave pairing
	Real p_p;				// p-wave pairing

	// The iteration variable.
	ImGreen gfloc;			// if in mode = 0, the APIzen iteration on G (default).
	ImGreen seloc;				// if in mode = 1, the APIzen iteration on self-energy.
private:
	bool converged() const {
		const Real dev = DEV(se_err);
		if(se_err[se_err.size()-2]<1.E-5 && se_err[se_err.size()-2]<se_err[se_err.size()-1]) return true;
		if(se_err[se_err.size()-1]<1.E-6) return true;
		if (dev > 1.E-4) { return false; }
		else if (dev > 1.E-10) {
			for_Int(i, 1, se_err.size()) {
				if (se_err[0] > se_err[i]) { return false; }
			}
			return true;
		}
		else{ return true; }
	}
	void append_se_err(Real err) {
		for_Int(i, 1, se_err.size()) {
			se_err[i - 1] = se_err[i];
		}
		se_err[se_err.size() - 1] = err;
	}
	void print_log(const Str& lbl, std::ostream& os = std::cout) const {
		using namespace std;
		os << iofmt("sci");
		os << setw(4) << iter_cnt
			<< "  " << setw(w_Real) << se_err[se_err.size() - 1]
			<< "  " << setw(w_Real) << p_n
			<< "  " << setw(w_Real) << m_m
			<< "  " << setw(w_Real) << d_p
			<< "  " << setw(w_Real) << p_p
			// << "  " << setw(w_Real) << mdl.epsd
			// << "  " << setw(w_Real) << mdl.epsp
			<< "  " << present()
			<< "  " << lbl << endl;
	}
	void log(const Str& lbl) {
		if (mm) {
			print_log(lbl);
			OFS ofs(p.ofx + ".log.txt", std::ios::app);
			print_log(lbl, ofs);
			ofs.close();
		}
	}
	// ImGreen find_hb(const ImGreen& g0eff_inv) const;
    // ImGreen rotate_hb(const ImGreen &hb) const;
    // void update_phy();
    // void save_the_backup(Bath &bth, NORG &solver, Int iter_cnt = 999);
	void shift_norg_restrain(){p.change_the_norg_restrain_and_div(p.stage2_restrain,  p.stage2_control);}
	void norg_up_date_h0_with_solve_mutable(const MatReal& h0_i, NORG& stage1_norg, NORG& stage2_norg){
		stage1_norg.up_date_h0_to_solve(h0_i);	stage2_norg.uormat = stage1_norg.uormat;
		stage2_norg.up_date_h0_to_solve(h0_i);	stage1_norg.uormat = stage2_norg.uormat;
	}
  public:
	APIzen(const MyMpi& mm_i, const Prmtr& prmtr_i, const Model& mdl_i);
};
