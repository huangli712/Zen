#pragma once

/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2023
*/

#include "operator.h"
#include "specs.h"
#include "prmtr.h"
#include "densitymatrix.h"
#include "nocspace.h"
#include "norg.h"

//The class for the Occler(occupancy controler)

class Occler{
	const MyMpi &mm;						// parameters
	Prmtr &p;							    // parameters
    VecInt &nparticals;                     // temporyly storge the partical number
    Real np_energy_a, np_energy_p, np_energy_m;

public:
	// VEC<MatReal> uormat;					// unitary_orbital_rotation_matrix

private:
    // void initial_uormat();
    Int if_ground_state();

    Str list_nppso(const VecInt& nppso_i) const{
		Str temp;
		for (const auto& i : nppso_i.truncate(0,nppso_i.size()/2)) temp += "-" + STR(i);
		return temp;
	}
public:
	Occler(const MyMpi& mm_i, Prmtr& prmtr_i);
    NORG find_ground_state_partical(const MatReal& h0_i);
};