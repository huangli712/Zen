#pragma once

/*
coded by Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2021-02-25
*/

#include "specs.h"
#include "prmtr.h"
#include "green.h"
#include "bath.h"
#include "model.h"

// a basis state (i.e., a configuration) is represented by an integer.
// the lower 16 bits represent a spin-down configuration
// the higher 16 bits represent a spin-up configuration

class Impurity {
private:
public:
	const MyMpi& mm;		// parameters
	const Prmtr& p;			// parameters
	const Bath& bth;
    // const Model &mdl;

    const Int ni;				//number of impurity
	const Int ns;				// number of sites,ns=ni+nb
	const Int nb;				// number of bath sites

	MatReal h0;					//hopping factors
	VecReal pos_imp;   			//position of imp site

private:
	//hopping  factors;when bath parameters is unusual,we just need to modify bath.hop() 
	void set_factor();
public:
	Impurity(const MyMpi& mm_i, const Prmtr& prmtr_i, const Bath& bth_i, const Str& file = empty_str);
	void find_g0(Green& g0) const;
	
	void find_hb(Green& hb) const;
	
	MatReal find_hop_for_test() const;

	void update();

	// void save() const;

	// void read(const Str& file) const;
};
