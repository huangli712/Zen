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
    // Real sub_energy[0], sub_energy[2], sub_energy[1];
	VecReal sub_energy;

public:
	// VEC<MatReal> uormat;					// unitary_orbital_rotation_matrix

private:
    // void initial_uormat();
    Int if_ground_state();
    // bool if_ground_state(Int counter, VecReal sub_energy);

	VEC<VecInt> list_all_posible_nppsos(const VecInt& nppso_i, const VecInt& or_deg) const ;

    Str nppso_str(const VecInt& nppso_i) const{
		Str temp;
		for_Int(i, 0, nppso_i.size()) if(i%2==0) temp += "-" + STR(nppso_i[i]);
		return temp;
	}

	VEC<VEC<Int>> cart_product(const VEC<VEC<Int>> &v) const
	{
		VEC<VEC<Int>> s = {{}};
		for (const auto &u : v)
		{
			VEC<VEC<Int>> r;
			for (const auto &x : s)
			{
				for (const auto y : u)
				{
					r.push_back(x);
					r.back().push_back(y);
				}
			}
			s = move(r);
		}
		return s;
	}



public:
	Occler(const MyMpi& mm_i, Prmtr& prmtr_i);
    NORG find_ground_state_partical(const MatReal& h0_i);
    NORG find_ground_state_partical(const MatReal& h0_i, const VecInt& or_deg);


};