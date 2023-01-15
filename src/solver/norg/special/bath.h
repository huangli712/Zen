#pragma once

/*
coded by Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2021-02-19
*/

#include "specs.h"
#include "prmtr.h"
#include "green.h"
#include "hyberr.h"

class Bath {
public:
	const MyMpi& mm;
	const Prmtr& p;			// parameters
	Int nb;					// number of bath sites
	ImGreen hb;				// hybridization function
	UURandomSFMT uur;
	VecReal ose;			// on-site energies for bath sites
	VecReal hop;			// hopping coefficients between impurity site and bath sites
	VEC<VecReal> vec_ose,vec_hop;
private:
	void regularize_ose_hop() {
		slctsort(ose, hop);
		// after a unitary transformation, hop can always be 
		// non-negative when there is only one impurity site
		hop = ABS(hop);
	}
	void init_ose_hop() {
		uur(ose);
		ose -= 0.5;
		ose *= 4.;
		uur(hop);
		hop -= 0.5;
		hop *= 4.;
		regularize_ose_hop();
	}
	std::tuple<Real, VecReal, Int> bath_fit_contest(const VecReal& a0);
	VecReal next_initial_fitting_parameters(const VecReal& a0, const Int& ntry_fine, Int& itry);
	void perturb(VecReal& a, Real amplitude, Real exponent) {
		for_Int(i, 0, a.size()) {
			Real ra = uur() - 0.5;
			Real re = uur() - 0.5;
			a[i] = a[i] * (1. + amplitude * ra) + SIGN(std::pow(10., (8 * ABS(re) - exponent)), a[i]);
		}
	}
public:
	Bath(const MyMpi& mm_i, const Prmtr& prmtr_i);
	void bath_fit(const ImGreen& hb_i, Int iter);
	MatReal find_hop() const;
	void write_ose_hop(Int iter_cnt) const {
		using namespace std;
		OFS ofs_app(p.ofx + ".ose.hop.txt", std::ios::app);
		ofs_app << iofmt("sci");
		ofs_app << setw(4) << iter_cnt;
		for_Int(i, 0, nb) { ofs_app << "\t" << setw(w_Real) << ose[i]; }
		for_Int(i, 0, nb) { ofs_app << "\t" << setw(w_Real) << hop[i]; }
		ofs_app << endl;
		ofs_app.close();

		OFS ofs(iox + "zic" + prefill0(iter_cnt, 3) + ".ose.hop.txt");
		ofs << iofmt("sci");
		for_Int(i, 0, nb) {
			ofs << setw(3) << i << "\t" << setw(w_Real) << ose[i] << "\t" << setw(w_Real) << hop[i] << endl;;
		}
		ofs.close();
	}
};