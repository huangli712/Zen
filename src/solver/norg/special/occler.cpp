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
    Int counter(0);
    while(1){
        if(mm) std::cout << "The " << ++counter << "-th NORG begin" << std::endl;	// norg counter
        NORG a(mm, p);
        a.up_date_h0_to_solve(h0_i);            np_energy_b = a.groune_lst;
        // VecInt np_m(p.npartical), np_p(p.npartical); np_m -= 1; np_p += 1;
        if(mm) std::cout << "The " << ++counter << "-th NORG begin" << std::endl;	// norg counter
        p.templet_control[1]--;     p.templet_control[p.ndiv-1]++;  p.after_modify_prmtr(); p.recalc_partical_number();
        NORG a_m(mm, p);
        a_m.up_date_h0_to_solve(h0_i);          np_energy_m = a_m.groune_lst;
        p.templet_control[1]++;     p.templet_control[p.ndiv-1]--;  p.after_modify_prmtr(); p.recalc_partical_number();
        if(mm) std::cout << "The " << ++counter << "-th NORG begin" << std::endl;	// norg counter
        p.templet_control[1]++;     p.templet_control[p.ndiv-1]--;  p.after_modify_prmtr(); p.recalc_partical_number();
        NORG a_p(mm, p);
        a_p.up_date_h0_to_solve(h0_i);          np_energy_p = a_p.groune_lst;
        p.templet_control[1]--;     p.templet_control[p.ndiv-1]++;  p.after_modify_prmtr(); p.recalc_partical_number();
        Int check = if_ground_state();
        if (check == 0) return a;
        // if (check == 1 && counter == 3) return a_p;
        // if (check == 2 && counter == 3) return a_m;
        if (check == 3) {
            if(mm) WRN("the occler is not converged.");
            return a;
        }
    }
    ERR("There some thing wrong in Occler.cpp!")
}










//-----------------------------------------private-----------------------------------------


// void Occler::initial_uormat(){
// 	for_Int(i, 0, p.nI2B.size()) {
// 		MatReal temp(dmat(p.nI2B[i], 1.));
// 		uormat.push_back(std::move(temp));
// 	}
// }

Int Occler::if_ground_state(){
    if(np_energy_b <= np_energy_m && np_energy_b <= np_energy_p) return 0;
    else
    {
        if (np_energy_m >= np_energy_b && np_energy_b > np_energy_p) {
            p.templet_control[1]++;     p.templet_control[p.ndiv-1]--;  p.after_modify_prmtr(); p.recalc_partical_number();
            return 1;
        }
        if ( np_energy_m < np_energy_b && np_energy_b <= np_energy_p) {
            p.templet_control[1]--;     p.templet_control[p.ndiv-1]++;  p.after_modify_prmtr(); p.recalc_partical_number();
            return 2;
        };
        if ( np_energy_m > np_energy_b && np_energy_b > np_energy_p)
            ERR("The form of energy is no longer dominated by the quadratic form, and the Occler class needs to be changed.");
    }
    return 3;
}