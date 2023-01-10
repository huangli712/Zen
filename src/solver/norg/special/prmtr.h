#pragma once

/*
coded by Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China) date 2020
*/

#include "specs.h"

class Prmtr {
private:
	const Int np;					// number of mpi processes
public:
	// model related, model parameters
	Str project;					// project name
	Int norbs;						// number of orbitals with spin.	
    MatReal c2u;    				//unitary transformation for c2 symmetry

	
	// square lattice, model parameters
	Real hubbU;						// Hubbard interaction U
	Real mu;						// chemical potential
	VecReal t;						// The  hopping for sites,t[0] is (on-site energy-mu)
	//ReGreen
	Real freq_upp;					// upper bound of real frequency, 1.5 * (band upper bound) suggested
	Real freq_low;					// lower bound of real frequency, 1.5 * (band lower bound) suggested
	Real dlt_freq;					// delta frequency, freq = n * dlt_freq, freq[i] = freq_low + dlt_freq * i;
	Real eta_freq;					// eta of omega + i * eta
	Int nfreq;						// number of real frequencies = 1 + Int_ROUND((freq_upp - freq_low) / dlt_freq)	// dmft related
    VecCmplx Re_z;		

	//ImGreen		
	Real unit_omg;					// unit imaginary frequency, omg_n = (2 n + 1) unit_omg, 0.01 or 0.02 suggested
	Real max_omg;					// imaginary frequency cutoff, 4 * (half bandwidth) suggested
	Int num_omg;					// number of positive imaginary frequencies
    VecCmplx Im_z;		

	//fitting related		
	Real fit_max_omg;				// imaginary frequency upper bound used for bath fitting
	Int fit_num_omg;				// number of positive imaginary frequencies used for bath fitting
	Real fit_pow;					// 1 / (omg_n + fit_rsd) ^ fit_pow is the weight used in bath fitting
	Real fit_rsd;					// 1 / (omg_n + fit_rsd) ^ fit_pow is the weight used in bath fitting

	//iter related
	bool imp_backup;				// imp backup
	Int iter_max;					// DMFT max number of iterations
	Int gauss_n_max;				// max gauss_n for Green function integration
	Int gauss_n_min;				// min gauss_n for Green function integration

	// New prmtr for NORG solver.
	Int norg_sets;					// number of norg sets (spinless orbital).
	Int iter_max_norg;				// the max NORG iteration times.
	mutable VecInt nI2B;			// number of bath sites to each impurity.
	mutable VecInt npartical;		// number of particals.
	mutable Int nbath;				// number of bath sites, must be an integer multiple of 4 
	mutable Int norbit;				// number of NORG orbital(imp + bath).

	mutable MatInt control_divs;	// to set the number of division and the shortcut restratin.
	mutable VecInt templet_restrain;// to set the restrain templet one for all distribute;
	mutable VecInt templet_control;	// to set the control templet one for all distribute;
	mutable VecInt stage2_restrain;	// to set the stage2 restrain one for all distribute;
	mutable VecInt stage2_control;	// to set the stage2 control one for all distribute;
	mutable Int ndiv;				// the the divsion's nubmer.

	
	mutable Str nooc_mode;

	
	// derived
	Str ofx;				// output filename prefix

private:
	void set_inert_values();
	void set_values();
	void derive();
	void check_consistency();
	void print(std::ostream &os, const Str &var, const Str &val, const Str &comment) const {
		using namespace std;
		Str true_var = var == "\"\"" ? "" : var;
		os << rght_justify(true_var, 16) << (true_var == "" ? "   " : " = ") << left_justify(val, w_Real)
			<< "    # " + comment << endl;
	}
public:
	Prmtr(const MyMpi& mm);
	void after_modify_prmtr() const;
	void print(std::ostream &os = std::cout) const;
	void change_the_norg_restrain_and_div(VecInt new_restrain, VecInt new_control) const;

	// Imag frequency
	Real Imomg(Int n) const { return imag(Im_z[n]); }
	Cmplx Imz(Int n) const { return Im_z[n]; }
	// Real frequency
	Real Reomg(Int n) const { return real(Re_z[n]); }
	Cmplx Rez(Int n) const { return Re_z[n]; }
};