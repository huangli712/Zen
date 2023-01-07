/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

#include "crrltfun.h"
using namespace std;

CrrltFun::CrrltFun(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& old_nosp_i, const NocSpace& main_nosp, const VecReal& vgs_i, const Int position)
    :Operator(mm_i, prmtr_i, main_nosp),
    old_nosp(old_nosp_i), new_nosp(main_nosp), crtorann(main_nosp.nspa - old_nosp.nspa),
    ex_state(project_uplwer_parical_space(vgs_i, crtorann, position))
{}

CrrltFun::CrrltFun(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& old_nosp_i, const NocSpace& main_nosp, const Tab& table, const VecReal& vgs_i, const Int position)
    :Operator(mm_i, prmtr_i, main_nosp, table),
    old_nosp(old_nosp_i), new_nosp(main_nosp), crtorann(main_nosp.nspa - old_nosp.nspa),
    ex_state(project_uplwer_parical_space(vgs_i, crtorann, position))
{}

CrrltFun::CrrltFun(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& old_nosp_i, const NocSpace& main_nosp, const Tab& table, const VecReal& vgs_i, const Int pos_in_set, const Int pos_in_div)
    :Operator(mm_i, prmtr_i, main_nosp, table),
    old_nosp(old_nosp_i), new_nosp(main_nosp), crtorann(main_nosp.nspa - old_nosp.nspa),
    ex_state(project_uplwer_parical_space(vgs_i, crtorann, pos_in_set, pos_in_div))
{}


// ?CrrltFun::CrrltFun(const NocSpace& old_nosp_i, const Operator& main_op, const VecReal& vgs_i, const Int position = 0)
// ?    :Operator(main_op),
// ?    old_nosp(old_nosp_i), new_nosp(main_op.scsp), crtorann(main_op.scsp.nspa - old_nosp.nspa),
// ?    ex_state(project_uplwer_parical_space(vgs_i, crtorann, position))
// ?{}

// (Deactivate) Only for test
ImGreen CrrltFun::find_density_density_correlation_function(const Real &ge0) {

    Real upper_fraction(DOT(ex_state, ex_state));
    Vectrdgnl trdgnl(find_trdgnl(ex_state, -1));
    // if (mm) WRN("Finish the find_trdgnl");
    VecReal zerod(get<0>(trdgnl).size(), 0.);
    VecReal zeros(get<1>(trdgnl).size(), 0.);
    VecCmplx a(cmplx(get<0>(trdgnl), zerod));
    VecCmplx b(cmplx(get<1>(trdgnl), zeros));
    Int n(a.size());
    ImGreen green(1, p);
    for_Int(w, 0, green.nomgs) {
        //Cmplx z(cmplx(ge0, green.omg(w)));
        Cmplx z = green.z(w) + ge0;
        Cmplx c(z + a[n - 2] - b[n - 2] * b[n - 2] / (z + a[n - 1]));
        for (Int i = n; i >= 3; i--) c = z + a[i - 3] - b[i - 3] * b[i - 3] / c;
        Cmplx under_fraction(c);
        Cmplx gaz(upper_fraction / under_fraction);
        green[w][0][0] = gaz;
    }
    return green;
}


void CrrltFun::find_gf_greater(const Real& ge0, Green &g0) 
{
    Real upper_fraction(DOT(ex_state, ex_state));    // if(mm)WRN(NAV(upper_fraction));
    //VECtrdgnl trdignl(find_trdgnl_first(ex_state));
    VEC<Real> ltd;	        // diagonal elements 
    VEC<Real> lt_sd;	    // sub-diagonal elements
    const SparseMatReal sep_h = find_hmlt(table);
    VecReal v0(ex_state);
    v0.normalize();
    VecReal v0_saved(v0);
    VecReal v1(sep_h * v0);
    Real a_i(0.), b_i(0.);
    a_i = DOT(v1, v0);
    ltd.push_back(a_i);
    if(g0.type_info() == STR("ImGreen")){
        for_Int(i, 0, 60) {
        find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
        ltd.push_back(a_i); lt_sd.push_back(b_i);
        }
    }
    if(g0.type_info() == STR("ReGreen")){
        for_Int(i, 0, 3000) {
        find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
        ltd.push_back(a_i); lt_sd.push_back(b_i);
        }
    }
    // ImGreen last_green(1, p);
    Cmplx zero(0.,0.);
    VecCmplx green_pre(g0.nomgs, zero);
    VecCmplx green_error(g0.nomgs, zero);
    while (true) {
        find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
        ltd.push_back(a_i); lt_sd.push_back(b_i);
        VecReal a(ltd), b(lt_sd);
        Int n(a.size());
        //ImGreen green_i(1, p);
        for_Int(w, 0, g0.nomgs) {
            Cmplx z = g0.z(w) + ge0;
            Cmplx c(z - a[n - 2] - b[n - 2] * b[n - 2] / (z - a[n - 1]));
            for (Int i = n; i >= 3; i--) c = z - a[i - 3] - b[i - 3] * b[i - 3] / c;
            //Cmplx under_fraction(c);
            Cmplx gaz = upper_fraction / c;
            green_error[w] = gaz - green_pre[w];
            green_pre[w] = gaz;
        }
        if (ABS(SUM(green_error)) < 1.E-10 * g0.nomgs) break;
    }
    for_Int(w, 0, g0.nomgs) g0[w][0][0] += green_pre[w];
    if (mm) PIO("The size of a and b in greaer:" + NAV2(ltd.size(), lt_sd.size()));
    // return last_green;
}

void CrrltFun::find_gf_lesser(const Real& ge0, Green &g0) 
{
    Real upper_fraction(DOT(ex_state, ex_state));
    //VECtrdgnl trdignl(find_trdgnl_first(ex_state));
    VEC<Real> ltd;	        // diagonal elements 
    VEC<Real> lt_sd;	    // sub-diagonal elements
    const SparseMatReal sep_h = find_hmlt(table);
    VecReal v0(ex_state);
    v0.normalize();
    VecReal v0_saved(v0);
    VecReal v1(sep_h * v0);
    Real a_i(0.), b_i(0.);
    a_i = DOT(v1, v0);
    ltd.push_back(a_i);
    if(g0.type_info() == STR("ImGreen")){
        for_Int(i, 0, 60) {
        find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
        ltd.push_back(a_i); lt_sd.push_back(b_i);
        }
    }
    if(g0.type_info() == STR("ReGreen")){
        for_Int(i, 0, 3000) {
        find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
        ltd.push_back(a_i); lt_sd.push_back(b_i);
        }
    }
    // ImGreen last_green(1, p);
    Cmplx zero(0.,0.);
    VecCmplx green_pre(g0.nomgs, zero);
    VecCmplx green_error(g0.nomgs, zero);
    while (true) {
        find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
        ltd.push_back(a_i); lt_sd.push_back(b_i);
        VecReal a(ltd), b(lt_sd);
        Int n(a.size());
        //ImGreen green_i(1, p);
        for_Int(w, 0, g0.nomgs) {
            Cmplx z = g0.z(w) - ge0;
            Cmplx c(z + a[n - 2] - b[n - 2] * b[n - 2] / (z + a[n - 1]));
            for (Int i = n; i >= 3; i--) c = z + a[i - 3] - b[i - 3] * b[i - 3] / c;
            //Cmplx under_fraction(c);
            Cmplx gaz = upper_fraction / c;
            green_error[w] = gaz - green_pre[w];
            green_pre[w] = gaz;
        }
        if (ABS(SUM(green_error)) < 1.E-10 * g0.nomgs) break;
    }
    for_Int(w, 0, g0.nomgs) g0[w][0][0] += green_pre[w];
    if (mm) PIO("The size of a and b in lesser:" + NAV2(ltd.size(), lt_sd.size()));
	// return last_green;
}


// void CrrltFun::find_gf_greater(const Real& ge0, Green &g0, Int kind) 
// {
//     Real upper_fraction(DOT(ex_state, ex_state));
//     //VECtrdgnl trdignl(find_trdgnl_first(ex_state));
//     VEC<Real> ltd;	        // diagonal elements 
//     VEC<Real> lt_sd;	    // sub-diagonal elements
//     const SparseMatReal sep_h = find_hmlt(table);
//     Cmplx zero(0.,0.);
//     VecCmplx green_pre(g0.nomgs, zero);
//     VecCmplx green_error(g0.nomgs, zero);
//     if(kind == 0){
//         VecReal v0(ex_state);
//         v0.normalize();
//         VecReal v0_saved(v0);
//         VecReal v1(sep_h * v0);
//         Real a_i(0.), b_i(0.);
//         a_i = DOT(v1, v0);
//         ltd.push_back(a_i);
//         for_Int(i, 0, 3000) {
//             find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
//             ltd.push_back(a_i); lt_sd.push_back(b_i);
//         }
//         // ImGreen last_green(1, p);
//         while (true) {
//             find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
//             ltd.push_back(a_i); lt_sd.push_back(b_i);
//             VecReal a(ltd), b(lt_sd);
//             Int n(a.size());
//             //ImGreen green_i(1, p);
//             for_Int(w, 0, g0.nomgs) {
//                 Cmplx z = g0.z(w, kind) + ge0;
//                 Cmplx c(z - a[n - 2] - b[n - 2] * b[n - 2] / (z - a[n - 1]));
//                 for (Int i = n; i >= 3; i--) c = z - a[i - 3] - b[i - 3] * b[i - 3] / c;
//                 //Cmplx under_fraction(c);
//                 Cmplx gaz = upper_fraction / c;
//                 green_error[w] = gaz - green_pre[w];
//                 green_pre[w] = gaz;
//             }
//             if (ABS(SUM(green_error)) < 1.E-10 * g0.nomgs) break;
//             a_saved.reset(a); b_saved.reset(b);
//         }
//     }
//     if(kind > 0){
//         VecReal a(a_saved), b(b_saved);
//         Int n(a.size());
//         for_Int(w, 0, g0.nomgs) {
//             Cmplx z = g0.z(w, kind) + ge0;
//             Cmplx c(z - a[n - 2] - b[n - 2] * b[n - 2] / (z - a[n - 1]));
//             for (Int i = n; i >= 3; i--) c = z - a[i - 3] - b[i - 3] * b[i - 3] / c;
//             //Cmplx under_fraction(c);
//             Cmplx gaz = upper_fraction / c;
//             green_pre[w] = gaz;
//         }
//     }
//     for_Int(w, 0, g0.nomgs) g0[w][0][0] += green_pre[w];
//     if (mm && kind == 0) PIO("The size of a and b in greaer:" + NAV2(a_saved.size(), b_saved.size()));
//     // return last_green;
// }

// void CrrltFun::find_gf_lesser(const Real& ge0, Green &g0, Int kind) 
// {
//     Real upper_fraction(DOT(ex_state, ex_state));
//     //VECtrdgnl trdignl(find_trdgnl_first(ex_state));
//     VEC<Real> ltd;	        // diagonal elements 
//     VEC<Real> lt_sd;	    // sub-diagonal elements
//     const SparseMatReal sep_h = find_hmlt(table);
//     Cmplx zero(0.,0.);
//     VecCmplx green_pre(g0.nomgs, zero);
//     VecCmplx green_error(g0.nomgs, zero);
//     if(kind == 0){
//         VecReal v0(ex_state);
//         v0.normalize();
//         VecReal v0_saved(v0);
//         VecReal v1(sep_h * v0);
//         Real a_i(0.), b_i(0.);
//         a_i = DOT(v1, v0);
//         ltd.push_back(a_i);
//         for_Int(i, 0, 3000) {
//             find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
//             ltd.push_back(a_i); lt_sd.push_back(b_i);
//         }
//         // ImGreen last_green(1, p);
//         while (true) {
//             find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
//             ltd.push_back(a_i); lt_sd.push_back(b_i);
//             VecReal a(ltd), b(lt_sd);
//             Int n(a.size());
//             //ImGreen green_i(1, p);
//             for_Int(w, 0, g0.nomgs) {
//                 Cmplx z = g0.z(w, kind) - ge0;
//                 Cmplx c(z + a[n - 2] - b[n - 2] * b[n - 2] / (z + a[n - 1]));
//                 for (Int i = n; i >= 3; i--) c = z + a[i - 3] - b[i - 3] * b[i - 3] / c;
//                 //Cmplx under_fraction(c);
//                 Cmplx gaz = upper_fraction / c;
//                 green_error[w] = gaz - green_pre[w];
//                 green_pre[w] = gaz;
//             }
//             if (ABS(SUM(green_error)) < 1.E-10 * g0.nomgs) break;
//             a_saved.reset(a); b_saved.reset(b);
//         }
//     }
//     if(kind > 0){
//         VecReal a(a_saved), b(b_saved);
//         Int n(a.size());
//         for_Int(w, 0, g0.nomgs) {
//             Cmplx z = g0.z(w, kind) - ge0;
//             Cmplx c(z + a[n - 2] - b[n - 2] * b[n - 2] / (z + a[n - 1]));
//             for (Int i = n; i >= 3; i--) c = z + a[i - 3] - b[i - 3] * b[i - 3] / c;
//             //Cmplx under_fraction(c);
//             Cmplx gaz = upper_fraction / c;
//             green_pre[w] = gaz;
//         }
//     }
//     for_Int(w, 0, g0.nomgs) g0[w][0][0] += green_pre[w];
//     if (mm && kind == 0) PIO("The size of a and b in lesser:" + NAV2(a_saved.size(), b_saved.size()));
//     // return last_green;
// }

//---------------------------------------------Private function---------------------------------------------

// (Deactivate) 
CrrltFun::Vectrdgnl CrrltFun::find_trdgnl(const VecReal& initial_vector, const Int crtann) {
    // for tridiagonal matrix
    VEC<Real> ltd;	    // diagonal elements 
    VEC<Real> lt_sd;	// sub-diagonal elements
    SparseMatReal sep_h = find_hmlt(table);
    VecReal v0(initial_vector);
    v0.normalize();
    VecReal v0_saved(v0);
    //DBG(NAV(v0_saved.truncate(0, 1000)));
    VecReal v1(sep_h * v0);
    Idx k = 0;
    Real a(0.), b(0.);
    while (true)
    {
        a = DOT(v1, v0);
        ltd.push_back(a);
        if (k >= 30) if (k > v0_saved.size() - 20) break;
        k++;
        v1 -= a * v0;
        v1 -= DOT(v1, v0_saved) * v0_saved;
        b = v1.norm();
        lt_sd.push_back(b);
        v1 *= (1 / b);
        SWAP(v0, v1);
        v1 *= -b;
        v1 += sep_h * v0;
    }
    VecReal va_i(ltd), vb_i(lt_sd);
#ifdef _ASSERTION_
    if (va_i.size() - 1 != vb_i.size())ERR("Wrong with count orbit idx: " + NAV2(va_i.size(), vb_i.size()));
#endif
    return std::make_tuple(va_i, vb_i);
}

void CrrltFun::find_trdgnl_one_step(const VecReal& initial_vector, VecReal& v0, VecReal& v1, Real&a, Real& b, const SparseMatReal& sep_h){
    v1 -= a * v0;
    b = v1.norm();
    v1 -= DOT(v1, initial_vector) * initial_vector;
    v1.normalize();
    SWAP(v0, v1);
    v1 *= -b;
    v1 += sep_h * v0;
    a = DOT(v1, v0);
}

VecReal CrrltFun::project_uplwer_parical_space(const VecReal &initial_vector, const Int crtann, const Int norg_set, const Int orbit_pos_in_div) const
// By using the C^+ or C on a spinless orbital(for the spinless reason ONE normal orbital has two spinless orbits).
// the norg_set only suppose for the even number.
{
    // clock_t t_find_Newstate; TIME_BGN("find_Newstate" + NAV(mm.id()), t_find_Newstate);
    VecReal ex_state_part(new_nosp.dim, 0.);
    VecPartition row_H(mm.np(), mm.id(), old_nosp.dim);
    for_Int(h_i, row_H.bgn(), row_H.end())
    {
        const Int subnosp(old_nosp.wherein_NocSpace(h_i));
        MatInt new_nospdiv = old_nosp.div[subnosp];
        if (crtann == -1) --new_nospdiv[norg_set][0];
        if (crtann == +1) ++new_nospdiv[norg_set][0];
        // if (new_nosp.ifin_NocSpace(new_nospdiv)) {
        if (new_nosp.ifin_NocSpace(new_nospdiv, new_nosp.nppso)) {
            const ComDivs group(h_i - old_nosp.idx_div[subnosp], (old_nosp.div[subnosp]), (old_nosp.sit_mat), true);
            VecOnb exd_cf = group.cf;
            if(if_in_this_orbital(exd_cf, crtann, norg_set, orbit_pos_in_div)) {
                if (crtann == -1) exd_cf[norg_set * new_nosp.ndivs] = exd_cf[norg_set * new_nosp.ndivs].ann(orbit_pos_in_div);
                if (crtann == +1) exd_cf[norg_set * new_nosp.ndivs] = exd_cf[norg_set * new_nosp.ndivs].crt(orbit_pos_in_div);
                const ComDivs b(exd_cf, new_nospdiv, old_nosp.sit_mat);
                Int begin_idx(-1);
                begin_idx = new_nosp.divs_to_idx.at(new_nospdiv.vec().string());
                if (begin_idx == -1) ERR("wrong with ex_state" + NAV2(group.ne, new_nospdiv));
                ex_state_part[begin_idx + b.idx] = exd_cf[norg_set * new_nosp.ndivs].sgn(orbit_pos_in_div) * initial_vector[h_i];
            }
        }
    }
    VecReal ex_state(mm.Allreduce(ex_state_part));
    // TIME_END("find_Newstate" + NAV(mm.id()), t_find_Newstate);
    return ex_state;
}

bool CrrltFun::if_in_this_orbital(const VecOnb &exd_cf, const Int crtann, const Int norg_set, const Int orbit_pos_in_div) const {
    if (crtann == -1) return exd_cf[norg_set * new_nosp.ndivs].isocc(orbit_pos_in_div);
    if (crtann == +1) return exd_cf[norg_set * new_nosp.ndivs].isuno(orbit_pos_in_div);
    ERR("Some thing wrong with this function, if_in_this_orbital ()!")
}

//
//
//VecReal CrrltFun::project_uplwer_parical_space(const VecReal& initial_vector, const Int crtann, const Int imp_div, const Int orbit_pos_in_div)const
//// By using the C^+ or C on a spinless orbital(for the spinless reason ONE normal orbital has two spinless orbits).
//// the imp_div only suppose for the even number.
//{
//    clock_t t_find_lowerhmlt;
//    TIME_BGN("find_lowerhmlt" + NAV(mm.id()), t_find_lowerhmlt);
//    VecPartition row_H(mm.np(), mm.id(), new_nosp.dim);
//    VecReal ex_state_part(row_H.len(), 0.);
//    for_Int(h_i, row_H.bgn(), row_H.end())
//    {
//        const Int subscsp(new_nosp.wherein_NocSpace(h_i));
//        const ComDivs cfg(h_i - new_nosp.idx_div[subscsp], Mat2Vec(new_nosp.div[subscsp]), Mat2Vec(new_nosp.sit_mat), true);
//        //spinless term
//        if (crtann == -1) {
//            //for_Int(updw, imp_div, imp_div + 2) {// use the now cfg to find the idx
//                if (cfg.cf[0 + imp_div * new_nosp.ndivs].isuno(orbit_pos_in_div)) {
//                    MatInt newdiv(new_nosp.div[subscsp]);
//                    ++newdiv[imp_div][0];
//                    if (old_nosp.ifin_NocSpace_for_green(newdiv)) {
//                        VecOnb find_old_cfg(cfg.cf);
//                        find_old_cfg[0 + imp_div * new_nosp.ndivs] = find_old_cfg[0 + imp_div * new_nosp.ndivs].crt(orbit_pos_in_div);
//                        Int begin_idx(-1);
//                        for_Int(j, 0, old_nosp.div.size()) if (old_nosp.div[j] == newdiv) { begin_idx = old_nosp.idx_div[j]; break; }
//                        if (begin_idx == -1)ERR("wrong with ex_state" + NAV2(cfg.ne, newdiv));
//                        ComDivs n(find_old_cfg, Mat2Vec(newdiv), Mat2Vec(old_nosp.sit_mat));
//                        ex_state_part[h_i - row_H.bgn()] += initial_vector[begin_idx + n.idx];
//                    }
//                }
//            //}
//        }
//        else if (crtann == +1) {
//            //for_Int(updw, imp_div, imp_div + 2) {// use the now cfg to find the idx
//                if (cfg.cf[0 + imp_div * new_nosp.ndivs].isocc(orbit_pos_in_div)) {
//                    MatInt newdiv(new_nosp.div[subscsp]);
//                    --newdiv[imp_div][0];
//                    if (old_nosp.ifin_NocSpace_for_green(newdiv)) {
//                        VecOnb find_old_cfg(cfg.cf);
//                        find_old_cfg[0 + imp_div * new_nosp.ndivs] = find_old_cfg[0 + imp_div * new_nosp.ndivs].ann(orbit_pos_in_div);
//                        Int begin_idx(-1);
//                        for_Int(j, 0, old_nosp.div.size()) if (old_nosp.div[j] == newdiv) { begin_idx = old_nosp.idx_div[j]; break; }
//                        if (begin_idx == -1)ERR("wrong with ex_state" + NAV2(cfg.ne, newdiv));
//                        ComDivs n(find_old_cfg, Mat2Vec(newdiv), Mat2Vec(old_nosp.sit_mat));
//                        ex_state_part[h_i - row_H.bgn()] += initial_vector[begin_idx + n.idx];
//                    }
//                }
//            //}
//        }
//    }
//    VecReal ex_state(new_nosp.dim, 0.);
//    ex_state = mm.Allgatherv(ex_state_part, row_H);
//    TIME_END("find_lowerhmlt" + NAV(mm.id()), t_find_lowerhmlt);
//    return ex_state;
//}