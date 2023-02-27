#include "impurity.h"

Impurity::Impurity(const MyMpi &mm_i, const Prmtr &prmtr_i, const Bath &bth_i, const Str& file)
    : mm(mm_i), p(prmtr_i), bth(bth_i), nb(p.nbath), ni(p.norbs), ns(p.norbit), pos_imp(p.norbs), h0(p.norbit, p.norbit, 0.)
{
    // if (!file.empty()) read(file);
    set_factor();
}

using namespace std;

// rely on imp model's frame
void Impurity::find_g0(Green &g0) const {

    MatCmplx Z(ns, ns);
    for_Int(n, 0, g0.nomgs) {
        Z = g0.z(n);
        MatCmplx g0_z = matinvlu(Z - cmplx(h0));
        for_Int(i, 0, p.nband) {
            g0[n][i][i] = g0_z[pos_imp[2 * i]][pos_imp[2 * i]];
        }
    }
}

// rely on imp model's frame
void Impurity::find_all_g0(Green &g0) const {

    MatCmplx Z(ns, ns);
    for_Int(n, 0, g0.nomgs) {
        Z = g0.z(n);
        MatCmplx g0_z = matinvlu(Z - cmplx(h0));
        for_Int(i, 0, p.norbit) {
            g0[n][i][i] = g0_z[i][i];
        }
    }
}

//rely on imp model's frame
void Impurity::find_hb(Green &hb) const {
    for_Int(i,0,p.nband){
        const Int nb_i = p.nI2B[2*i];
        for_Int(n,0,hb.nomgs){
            const VecCmplx Z(nb_i, hb.z(n));
            VecCmplx S = INV(Z - cmplx(bth.vec_ose[i]));
            VecCmplx V = cmplx(bth.vec_hop[i]);
            hb[n][i][i] = -SUM(V * S * V.co());
            //hb[n][i][i] = DOT(cmplx(bth.vec_hop[i]), S * cmplx(bth.vec_hop[i]));
        }
    }
}

MatReal Impurity::find_hop_for_test() const
{
    VecReal hop({ -0.312546, 0.159346, -0.0619612, -0.312546, 0.159346 });
    VecReal ose({ -0.474855,0.0667285, 0,           0.474855, -0.0667285 });
    MatReal h0(0, 0, 0.);
    for_Int(i, 0, p.nband) {
        const Int nb_i = p.nI2B[i * 2];
        MatReal h0_i(1 + nb_i, 1 + nb_i, 0.);
        for_Int(j, 0, nb_i) {
            h0_i[0][j + 1] = hop[j];
            h0_i[j + 1][0] = hop[j];
            h0_i[j + 1][j + 1] = ose[j];
        }
        h0_i.reset(direct_sum(h0_i, h0_i));
        h0.reset(direct_sum(h0, h0_i));
    }
    return h0;
}

void Impurity::update() {
    set_factor();
}

void Impurity::write_H0info(const Bath &b, Int ndeg) const {
    OFS ofs;    ofs.open("h0.txt");
    using namespace std;
    if(ndeg > 0) for_Int(i, 0, ndeg)	{
        ofs << "The impurity for "<< i+1 << "-th degeneracy: " << "nmin: " << b.info[i][0] << " err: " << b.info[i][1] << " err_crv: " << b.info[i][2] << " err_reg: " << b.info[i][3] << " norm: " << b.info[i][4]<< "  " << endl;
    }
    else for_Int(i, 0, p.nband) {
        ofs << "The impurity for "<< i+1 << "-th band: " << "nmin: " << b.info[i][0] << " err: " << b.info[i][1] << " err_crv: " << b.info[i][2] << " err_reg: " << b.info[i][3] << " norm: " << b.info[i][4]<< "  " << endl;
    }
    for_Int(i, 0, p.nband)	{
        ofs << i+1 << "-th band: " << endl;
        Int begin(i*2 * (p.nI2B[i*2] + 1)), end((i*2 + 1) * (p.nI2B[i*2] + 1));
        ofs << iofmt() << h0.truncate(begin, begin, end, end) << endl;
    }
}

//---------------------------------------------Private function---------------------------------------------


void Impurity::set_factor() {
    // set hyb part and bath part
    h0 = bth.find_hop();
    
    // // h0 = find_hop_for_test();
    // set imp part
    MatReal h0loc(ni,ni,0.);
    for_Int(i, 0, ni) h0loc[i][i] = p.eimp[i] - p.mu;
    
    // find i_th imp in which site
    Int site = 0;
    for_Int(i, 0, ni) {
        pos_imp[i] = site;
        site += p.nI2B[i] + 1;
    }

    // set h0loc
    for_Int(i, 0, ni) {
        for_Int(j, 0, ni) {
            h0[pos_imp[i]][pos_imp[j]] = h0loc[i][j];
        }
    }
}