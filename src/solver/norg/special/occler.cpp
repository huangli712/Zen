/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

#include "occler.h"
using namespace std;

Occler::Occler(const MyMpi& mm_i, Prmtr& prmtr_i):
    mm(mm_i), p(prmtr_i), nparticals(p.npartical),
    np_energy_b(1.), np_energy_p(0.), np_energy_m(0.)
{
    // initial_uormat();
}



NORG Occler::find_ground_state_partical(const MatReal& h0_i){
    while(1){
        VecInt np_m(p.npartical), np_p(p.npartical); np_m -= 1; np_p += 1;
        if(mm) WRN(NAV3(p.npartical, np_m, np_p))
        NORG a(mm, p);
        a.up_date_h0_to_solve(h0_i);            np_energy_b = a.groune_lst;
        NORG a_p(mm, p, np_p);
        a_p.up_date_h0_to_solve(h0_i);          np_energy_p = a_p.groune_lst;
        NORG a_m(mm, p, np_m);
        a_m.up_date_h0_to_solve(h0_i);          np_energy_m = a_m.groune_lst;
        if (if_ground_state()) return a;
        // return a;
    }
}










//-----------------------------------------private-----------------------------------------


// void Occler::initial_uormat(){
// 	for_Int(i, 0, p.nI2B.size()) {
// 		MatReal temp(dmat(p.nI2B[i], 1.));
// 		uormat.push_back(std::move(temp));
// 	}
// }

bool Occler::if_ground_state(){
    if(np_energy_b <= np_energy_m && np_energy_b <= np_energy_p) return true;
    else
    {
        if (np_energy_b <= np_energy_m && np_energy_b > np_energy_p) nparticals += 1;
        if (np_energy_b > np_energy_m && np_energy_b <= np_energy_p) nparticals -= 1;
        if (np_energy_b > np_energy_m && np_energy_b > np_energy_p)
            ERR("The form of energy is no longer dominated by the quadratic form, and the Occler class needs to be changed.");
    }
    return false;
}