/*
coded by Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 20201028
*/

#include "model.h"

Model::Model(const MyMpi& mm_i, const Prmtr& prmtr_i) :
	mm(mm_i), p(prmtr_i), epsd(0.), epsp(0.)
{}

MatCmplx Model::hmlt0k(Real kx, Real ky) const
{
    MatCmplx tk(p.norbs, p.norbs, 0.);
	for_Int(i, 0, p.norbs) {
        tk[i][i] = -2 * p.t[3] * (cos(kx) + cos(ky)) + p.t[0];
    }
	tk[0][1] = -p.t[1] * (1. + exp(-I * kx));
    tk[1][0] = -p.t[1]* (1. + exp(I * kx));
    tk[0][2] = -p.t[1]* (1. + exp(-I * ky));
    tk[2][0] = -p.t[1]* (1. + exp(I * ky));
    tk[1][3] = -p.t[1]* (1. + exp(-I * ky));
    tk[3][1] = -p.t[1]* (1. + exp(I * ky));
    tk[2][3] = -p.t[1]* (1. + exp(-I * kx));
    tk[3][2] = -p.t[1]* (1. + exp(I * kx));

    // tk[0][1] = -p.t[1] * 2 * cos(kx);
    // tk[1][0] = -p.t[1] * 2 * cos(kx);
    // tk[0][2] = -p.t[1] * 2 * cos(ky);
    // tk[2][0] = -p.t[1] * 2 * cos(ky);
    // tk[1][3] = -p.t[1] * 2 * cos(ky);
    // tk[3][1] = -p.t[1] * 2 * cos(ky);
    // tk[2][3] = -p.t[1] * 2 * cos(kx);
    // tk[3][2] = -p.t[1] * 2 * cos(kx);

    tk[0][3] = -p.t[2] * (1. + exp(-I * kx)) * (1. + exp(-I * ky));
    tk[1][2] = -p.t[2] * (1. + exp(I * kx)) * (1. + exp(-I * ky));
    tk[2][1] = -p.t[2] * (1. + exp(-I * kx)) * (1. + exp(I * ky));
    tk[3][0] = -p.t[2] * (1. + exp(I * kx)) * (1. + exp(I * ky));
    return direct_sum(tk, -tk);
}

MatReal Model::hmlt0loc() const
{
    MatReal tloc(p.norbs, p.norbs, 0.);
    for_Int(i, 0, p.norbs) {
        tloc[i][i] = p.t[0];
        tloc[i][p.norbs - 1 - i] = -p.t[2];
    }
    tloc[0][1] = -p.t[1];
    tloc[1][0] = -p.t[1];
    tloc[0][2] = -p.t[1];
    tloc[2][0] = -p.t[1];
    tloc[1][3] = -p.t[1];
    tloc[3][1] = -p.t[1];
    tloc[2][3] = -p.t[1];
    tloc[3][2] = -p.t[1];
    MatReal h0loc = direct_sum(tloc, -tloc);
    return h0loc;
}

void Model::find_gfloc(const Green& seloc, Green & gfloc) const
{
	for_Int(n, 0, gfloc.nomgs) {
		const Cmplx z = gfloc.z(n);
		const Int gauss_n = Int_ROUND(p.gauss_n_min + (p.gauss_n_max - p.gauss_n_min) / (1 + 0.3 * n));
		VecPartition vp(mm.np(), mm.id(), gauss_n * gauss_n);
		MatCmplx mpi_g0loc = find_gfloc_part(vp.bgn(), vp.end(), gauss_n, z, seloc[n]);
		gfloc[n] = mm.Allreduce(mpi_g0loc);
		if (mm) {
//			std::cout << STR("find_g0loc, ") + NAV3(n, p.num_omg, gauss_n) << std::endl;
		}
	}
}


MatCmplx Model::find_gfloc_part(Int idx_bgn, Int idx_end, Int gauss_n, Cmplx z, const MatCmplx& seloc) const
{
    MatCmplx sum(2 * p.norbs, 2 * p.norbs, 0.);
    Gauss gauss(gauss_n, 0, -pi_Real, pi_Real);
    for_Int(idx, idx_bgn, idx_end) {
		Int j = idx / gauss_n;
		Int i = idx % gauss_n;
		Real ky = gauss.x[j];
		Real kx = gauss.x[i];
		MatCmplx G_k_omg = find_G_k_omg(kx, ky, z, seloc);
		sum += G_k_omg * Cmplx(gauss.w[i] * gauss.w[j]);
	}
	return sum * Cmplx(INV(4 * pi_Real * pi_Real));
}

MatCmplx Model::find_G_k_omg(Real kx, Real ky, Cmplx z, const MatCmplx& seloc) const {
    const MatCmplx h0k = hmlt0k(kx, ky);
    MatCmplx G_k_omg_inv(2 * p.norbs, 2 * p.norbs);
    G_k_omg_inv = dmat(2 * p.norbs, z) - h0k - seloc;
    MatCmplx G_k_omg = matinvlu(G_k_omg_inv);
	return G_k_omg;
}

void Model::find_gk(Real kx, Real ky, const Green &seloc, Green &gk) const {
    for_Int(n,0,seloc.nomgs){
        gk[n] = find_G_k_omg(kx, ky, seloc.z(n), seloc[n]);
    }
}