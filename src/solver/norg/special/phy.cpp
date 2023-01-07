#include "phy.h"


using namespace std;

Phy::Phy(const MyMpi& mm_i, const Prmtr& prmtr_i, const Model& mdl_i) :
	mm(mm_i), p(prmtr_i), mdl(mdl_i)
{
	{ mm.barrier(); SLEEP(1); }
   	Real M = find_M();  if(mm) std::cout << "M = " << M << endl;
    Cmplx D = find_D(); if(mm) std::cout << "D = " << D << endl;
    Cmplx P = find_P(); if(mm) std::cout << "P = " << P << endl;
}

Real Phy::find_M(){
    MatReal nloc_up = dmat(p.norbs,1.); 
    MatReal mloc_up = dmat(p.norbs, 1.);
    mloc_up[1][1] *= -1;
    mloc_up[2][2] *= -1;
    MatReal nloc = direct_sum(nloc_up, -nloc_up);
    MatReal mloc = direct_sum(mloc_up, mloc_up);

    ImGreen gfloc(p.norbs * 2, p, iox + "output_" + "gfloc" + ".txt");
    MatReal average = gfloc.particle_number();
    if(mm)   std::cout << "number of partical matrix is(n_down need 1-) : " << average << endl;
    Real m(0.);
    Real n(4.);
    for_Int(i, 0, p.norbs * 2) {
     	n += nloc[i][i] * average[i][i];
        m += mloc[i][i] * average[i][i];
    }
    if(mm)   std::cout << "number of partical : " << n / 4. << endl;
    return m / 4.;
}

Cmplx Phy::find_D() {
    ImGreen seloc(p.norbs * 2, p, iox + "output_" + "seloc" + ".txt");
    const Int gauss_n = Int_ROUND(p.gauss_n_min + (p.gauss_n_max - p.gauss_n_min) / (1 + 0.3));
    VecPartition vp(mm.np(), mm.id(), gauss_n * gauss_n);
    ImGreen gf(seloc);
    Gauss gauss(gauss_n, 0, -pi_Real, pi_Real);
    for_Int(n, 0, seloc.nomgs) {
        const Cmplx z = seloc.z(n);
        MatCmplx sum(seloc.norbs, seloc.norbs, 0.);
        for_Int(idx, vp.bgn(), vp.end()) {
            const Int j = idx / gauss_n;
            const Int i = idx % gauss_n;
            const Real ky = gauss.x[j];
            const Real kx = gauss.x[i];
            const MatCmplx dk = find_dk(kx, ky);
            MatCmplx gkomg = mdl.find_G_k_omg(kx, ky, z, seloc[n]);
            MatCmplx g_dk(gkomg);
            for_Int(row, 0, seloc.norbs) {
                for_Int(col,0,seloc.norbs){
                    g_dk[row][col] = gkomg[row][col] * dk[col][row];
                }
            }
            sum += g_dk * Cmplx(gauss.w[i] * gauss.w[j]);
            // if(mm) WRN(NAV(sum))
        }
        gf[n] = mm.Allreduce(sum * Cmplx(INV(4 * pi_Real * pi_Real)));
    }
    gf.write("gfd");
    MatReal npd = gf.particle_number() - dmat(p.norbs * 2, 0.5);
    if (mm) std::cout << "d-sc matrix is: " << npd << endl;
    return SUM(npd);
}

MatCmplx Phy::find_dk(Real kx, Real ky) {    
    MatCmplx d_k(p.norbs, p.norbs, 0.);
    d_k[0][1] = -(1. + exp(-I * kx));
    d_k[1][0] = -(1. + exp(I * kx));
    d_k[0][2] = (1. + exp(-I * ky));
    d_k[2][0] = (1. + exp(I * ky));
    d_k[1][3] = (1. + exp(-I * ky));
    d_k[3][1] = (1. + exp(I * ky));
    d_k[2][3] = -(1. + exp(-I * kx));
    d_k[3][2] = -(1. + exp(I * kx));

    MatCmplx dk(2 * p.norbs, 2 * p.norbs, 0.);
    for_Int(i, 0, p.norbs) {
        for_Int(j, 0, p.norbs) {
            dk[i][j + p.norbs] = d_k[i][j];
            dk[j + p.norbs][i] = d_k.ct()[j][i];
        }
    }
    return dk;
}

Cmplx Phy::find_P() {
    ImGreen seloc(p.norbs * 2, p, iox + "output_" + "seloc" + ".txt");
    const Int gauss_n = Int_ROUND(p.gauss_n_min + (p.gauss_n_max - p.gauss_n_min) / (1 + 0.3));
    VecPartition vp(mm.np(), mm.id(), gauss_n * gauss_n);
    ImGreen gf(seloc);
    Gauss gauss(gauss_n, 0, -pi_Real, pi_Real);
    for_Int(n, 0, seloc.nomgs) {
        const Cmplx z = seloc.z(n);
        MatCmplx sum(seloc.norbs, seloc.norbs, 0.);
        for_Int(idx, vp.bgn(), vp.end()) {
            const Int j = idx / gauss_n;
            const Int i = idx % gauss_n;
            const Real ky = gauss.x[j];
            const Real kx = gauss.x[i];
            const MatCmplx pk = find_pk(kx, ky);
            MatCmplx gkomg = mdl.find_G_k_omg(kx, ky, z, seloc[n]);
            MatCmplx g_pk(gkomg);
            for_Int(row, 0, seloc.norbs) {
                for_Int(col,0,seloc.norbs){
                    g_pk[row][col] = gkomg[row][col] * pk[col][row];
                }
            }
            sum += g_pk * Cmplx(gauss.w[i] * gauss.w[j]);
            // if(mm) WRN(NAV(sum))
        }
        gf[n] = mm.Allreduce(sum * Cmplx(INV(4 * pi_Real * pi_Real)));
    }
    gf.write("gfp");
    MatReal npp = gf.particle_number() - dmat(p.norbs * 2, 0.5);
    if (mm) std::cout << "pi matrix is: " << npp << endl;
    return SUM(npp);
}

MatCmplx Phy::find_pk(Real kx, Real ky) {    
    MatCmplx pi_k(p.norbs, p.norbs, 0.);
    pi_k[0][1] = -(1. + exp(-I * kx));
    pi_k[1][0] = (1. + exp(I * kx));
    pi_k[0][2] = (1. + exp(-I * ky));
    pi_k[2][0] = -(1. + exp(I * ky));
    pi_k[1][3] = -(1. + exp(-I * ky));
    pi_k[3][1] = (1. + exp(I * ky));
    pi_k[2][3] = (1. + exp(-I * kx));
    pi_k[3][2] = -(1. + exp(I * kx));

    MatCmplx pk(2 * p.norbs, 2 * p.norbs, 0.);
    for_Int(i, 0, p.norbs) {
        for_Int(j, 0, p.norbs) {
            pk[i][j + p.norbs] = pi_k[i][j];
            pk[j + p.norbs][i] = pi_k.ct()[j][i];
        }
    }
    return pk;
}

/*
Cmplx Phy::find_P() {
    const Int gauss_n = Int_ROUND(p.gauss_n_min + (p.gauss_n_max - p.gauss_n_min) / (1 + 0.3));
    VecPartition vp(mm.np(), mm.id(), gauss_n * gauss_n);
    Cmplx sum(0.);
    Gauss gauss(gauss_n, 0, -pi_Real, pi_Real);
    for_Int(idx, vp.bgn(), vp.end()) {
        Int j = idx / gauss_n;
        Int i = idx % gauss_n;
        Real ky = gauss.x[j];
        Real kx = gauss.x[i];
        Cmplx pk = find_pk(kx, ky);
        sum += pk * Cmplx(gauss.w[i] * gauss.w[j]);
    }
    Cmplx P = mm.Allreduce(sum);
    return P;
}



Cmplx Phy::find_pk(Real kx, Real ky) {
    MatCmplx pi_k(p.norbs, p.norbs, 0.);
    pi_k[0][1] = (1. - exp(-I * kx));
    pi_k[1][0] = (1. - exp(I * kx));
    pi_k[0][2] = -(1. - exp(-I * ky));
    pi_k[2][0] = -(1. - exp(I * ky));
    pi_k[1][3] = (1. - exp(-I * ky));
    pi_k[3][1] = (1. - exp(I * ky));
    pi_k[2][3] = -(1. - exp(-I * kx));
    pi_k[3][2] = -(1. - exp(I * kx));

    MatCmplx pik(2 * p.norbs, 2 * p.norbs, 0.);
    for_Int(i, 0, p.norbs) {
        for_Int(j, 0, p.norbs) {
            pik[i][j + p.norbs] = pi_k[i][j];
            pik[j + p.norbs][i] = pi_k.ct()[j][i];
        }
    }

    ImGreen seloc(p.norbs * 2, p, iox + "output_" + "seloc" + ".txt");
    ImGreen gk(p.norbs * 2, p);
    mdl.find_gk(kx, ky, seloc, gk);
    MatReal average = gk.particle_number();
    Cmplx pk(0.);
    for_Int(i, 0, p.norbs * 2) {
        for_Int(j, 0, p.norbs * 2) {
            pk += pik[i][j] * average[i][j];
        }
    }
    return pk;
}*/
