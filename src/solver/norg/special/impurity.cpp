#include "impurity.h"

Impurity::Impurity(const MyMpi &mm_i, const Prmtr &prmtr_i, const Bath &bth_i, const Model &mdl_i, const Str &file)
    : mm(mm_i), p(prmtr_i), bth(bth_i), mdl(mdl_i), nb(p.nbath), ni(2 * p.norbs), ns(2 * p.norbs + p.nbath),
      V(2 * p.norbs, p.nbath, 0.), E(p.nbath, p.nbath, 0.), h0(2 * p.norbs + p.nbath, 2 * p.norbs + p.nbath, 0.) {
    // find_ve();
    // find_ve_for_test();
    set_factor();
    // if(mm) WRN(NAV(h0))
}

using namespace std;

// rely on imp model's frame
void Impurity::find_g0(Green &g0) const {

    MatCmplx Z(ns, ns);
    // for_Int(kind, 0, p.a.size()) {
    for_Int(n, 0, g0.nomgs) {
        Z = g0.z(n);
        MatCmplx g0_z = matinvlu(Z - cmplx(h0));
        for_Int(i, 0, ni) {
            for_Int(j, 0, ni) {
                g0[n][i][j] = g0_z[i][j];
            }
        }
    }
}

//rely on imp model's frame
void Impurity::find_hb(Green &hb) const {
    for_Int(n, 0, hb.nomgs) {
        const MatCmplx Z=dmat(ni,hb.z(n));
        MatCmplx S = matinvlu(Z - cmplx(E));
        hb[n] = cmplx(V) * S * cmplx(V).co();
    }
}

// void Impurity::find_ve_for_test()
// {
//     MatReal c2vu(4, 4, 0.5);
//     c2vu[2][1] *= -1;
//     c2vu[3][1] *= -1;
//     c2vu[1][2] *= -1;
//     c2vu[3][2] *= -1;
//    	c2vu[1][3] *= -1;
//     c2vu[2][3] *= -1;
//     //nb=8*theta.size();
//     VecReal theta({-0.312546, 0.159346});
//     VecReal delta({-0.474855, 0.0667285});
//     VecReal epson({-0.474855, 0.0667285});
//     Int halfnb_qrt = theta.size();
//     VecReal ksi = SQRT(epson * epson + delta * delta);
//     VecReal u = delta * INV(SQRT(2. * ksi * (ksi - epson)));
//     VecReal v = SQRT(VecReal(halfnb_qrt, 1.) - u * u);
//     const MatReal unitm = dmat(4, 1.);
//     VEC<Real> E_up;
//     for_Int(n, 0, halfnb_qrt) {
//         MatReal Theta_n = theta[n] * unitm * c2vu;
//         VecReal un_vec({1., u[n], u[n], 1.});
//         VecReal vn_vec({0., -v[n], v[n], 0.});
//         VecReal en_vec({epson[n], ksi[n], ksi[n], epson[n]});
//         MatReal Theta = Theta_n * dmat(un_vec);
//         MatReal Delta = Theta_n * dmat(vn_vec);
//         for_Int(i, 0, 4) {
//             E_up.push_back(en_vec[i]);
//             for_Int(j, 0, 4) {
//                 Int col = j + 4 * n;
//                 V[i][col] = Theta[i][j];
//                 V[i + ni / 2][col + nb / 2] = -Theta[i][j];
//                 V[i][col + nb / 2] = -Delta[i][j];
//                 V[i + ni / 2][col] = -Delta[i][j];
//             }
//         }
//     }
//     E=direct_sum(dmat(VecReal(E_up)),-dmat(VecReal(E_up)));
// }

//---------------------------------------------Private function---------------------------------------------

//rely on the format of bath and model
void Impurity::set_factor() {
    // set hyb part and bath part
    // set imp part
    const MatReal h0loc = mdl.hmlt0loc();
    for_Int(i, 0, ni) {
        for_Int(j, 0, ni) {
            h0[i][j] = h0loc[i][j];
        }
    }

    for_Int(j, 0, nb) {
        const Int b = j + ni;
        for_Int(i, 0, ni) {
            h0[i][b] = V[i][j];
            h0[b][i] = V.ct()[j][i];
        }
        h0[b][b] = E[j][j];
    }
}


// void Impurity::find_ve(){
//     for_Int(i, 0, ni / 2) {
//         for_Int(j, 0, nb / 2) {
//             V[i][j] = bth.fvb[0][i][j] + bth.fvb[1][i][j];
//             V[i][j + nb / 2] = -bth.fvb[2][i][j] + bth.fvb[3][i][j];
//             V[i + ni / 2][j] = -bth.fvb[2][i][j] - bth.fvb[3][i][j];
//             V[i + ni / 2][j + nb / 2] = -bth.fvb[0][i][j] + bth.fvb[1][i][j];
//         }
//     }
//     V = direct_sum(p.c2u, p.c2u) * V;
//     E = direct_sum(bth.fvb[4] + bth.fvb[5], -bth.fvb[4] + bth.fvb[5]);
//     sort_v_and_e();
// }

// void Impurity::sort_v_and_e(){
//     VecReal e_temp(E.diagonal());
//     MatReal v_temp(V.ct());
//     slctsort(e_temp, v_temp);
//     E = dmat(e_temp); V = v_temp.ct();
// }