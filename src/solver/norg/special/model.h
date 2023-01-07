#pragma once

/*
coded by Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 20201011
*/

#include "specs.h"
#include "prmtr.h"
#include "green.h"

// dpsc model

class Model {
public:
	const MyMpi &mm;	// parameters
	const Prmtr &p;		// parameters
	Real epsd;			// current value for on-site energy for d orbital
	Real epsp;			// current value for on-site energy for p orbital
public:
	Model(const MyMpi& mm_i, const Prmtr &prmtr_i);
	MatCmplx hmlt0k(Real kx, Real ky) const;
	// h0 within the impurity cluster
	MatReal hmlt0loc() const;
    void find_gk(Real kx, Real ky, const Green &seloc, Green &gk) const;
    void find_gfloc(const Green &seloc, Green &gfloc) const;
    MatCmplx find_gfloc_part(Int idx_bgn, Int idx_end, Int gauss_n, Cmplx z, const MatCmplx &seloc) const;
	MatCmplx find_G_k_omg(Real kx, Real ky, Cmplx z, const MatCmplx& seloc) const;
};
