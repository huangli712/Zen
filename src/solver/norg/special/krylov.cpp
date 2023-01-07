/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

#include "krylov.h"
using namespace std;


Krylov::Krylov(const DensityMat& old_one_emat_i, const Int position_i) :DensityMat(old_one_emat_i),
ground_state_energy(old_one_emat_i.groundstate_energy), ground_state_i(old_one_emat_i.ground_state),
position(position_i)
{

}

ImGreen Krylov::find_gf_from_krylov_space(const Int crtann) {
    //Real upper_fraction(DOT(ex_state, ex_state));
    //VECtrdgnl trdignl(find_trdgnl_first(ex_state));
    VEC<Real> ltd;	        // diagonal elements 
    VEC<Real> lt_sd;	    // sub-diagonal elements
    VEC<Real> inner_m;	    // for the inner product of Krylov basis and ex_state.
    const SparseMatReal sep_h = find_hmlt(table);
    VecReal ex_state(sn_prtcl_ex_state(position, ground_state_i, crtann));
    //WRN(NAV3(ground_state_i, ex_state, crtann))
    VecReal v0(ex_state);
    v0.normalize();
    VecReal v0_saved(v0), v1(sep_h * v0);
    Real a_i(0.), b_i(0.);
    a_i = DOT(v1, v0);
    ltd.push_back(a_i); inner_m.push_back(DOT(v0, ex_state));
    //inner_m.push_back(DOT((v1*INV(v1.norm())), ex_state));
    for_Int(i, 0, 148) {
        find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
        lt_sd.push_back(b_i); ltd.push_back(a_i); inner_m.push_back(DOT(v0, ex_state));
    }
    Int counter(0);
    ImGreen last_green(1, p);
    while (true) {
        find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
        lt_sd.push_back(b_i); ltd.push_back(a_i), inner_m.push_back(DOT(v0, ex_state));
        VecReal m(inner_m), va(ltd), vb(lt_sd);
        VecCmplx green_error(last_green.nomgs);
        MatReal ev(ltd.size(), ltd.size(), 0.);
        Int info = trd_heevr_vd(va, vb, ev); ev = ev.tr();
        if (info != 0)ERR("the i-th parameter had an illegal value." + NAV3(info, va, vb));
        m = ev * m; 
        m.normalize();
        WRN(NAV(va[0]));
        if (crtann == +1)for_Int(w, 0, last_green.nomgs) {
            Cmplx gaz(0., 0.), z = last_green.z(w) + ground_state_energy;
            //for_Int(i, 0, va.size()) gaz += DOT(ev[i], m) * INV((z - va[i])) * DOT(m, ev[i]);
            for_Int(i, 0, va.size()) gaz += m[i] * INV((z - va[i])) * m[i];
            green_error[w] = gaz - last_green[w][0][0];
            last_green[w][0][0] = gaz;
        }
        if (crtann == -1)for_Int(w, 0, last_green.nomgs) {
            Cmplx gaz(0., 0.), z = last_green.z(w) - ground_state_energy;
            //for_Int(i, 0, va.size()) gaz += DOT(ev[i], m) * INV((z + va[i])) * DOT(m, ev[i]);
            for_Int(i, 0, va.size()) gaz += m[i] * INV((z + va[i])) * m[i];
            green_error[w] = gaz - last_green[w][0][0];
            last_green[w][0][0] = gaz;
        }
        WRN("The size of a and b in here:" + NAV2(va.size(), vb.size()));
        break;//TEST!!;
        
        counter++;
    }
    return last_green;
}

//---------------------------------------------Private function---------------------------------------------


void Krylov::find_trdgnl_one_step(const VecReal& initial_vector, VecReal& v0, VecReal& v1, Real& a, Real& b, const SparseMatReal& sep_h) {
    v1 -= a * v0;
    //v1 -= DOT(v1, initial_vector) * initial_vector;         // ? is necessary?
    b = v1.norm();
    v1 *= (1 / b);
    SWAP(v0, v1);
    v1 *= -b;
    v1 += sep_h * v0;
    a = DOT(v1, v0);
}