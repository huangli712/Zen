/*
coded by Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)date 2020-10-11
*/

#include "prmtr.h"

Prmtr::Prmtr(const MyMpi &mm) : np(mm.np())
{
    set_inert_values();
    set_values();
    derive();
    if (mm) print();
    if (mm) { OFS ofs(iox + "log.parameters.txt", std::ios::app);  print(ofs); }
    if (mm) { OFS ofss(cdr + "log.parameters.txt", std::ios::app);  print(ofss); }
}

void Prmtr::set_inert_values()
{
    project = "3band-SrVO3";
    nband = 3;                 
    norbs = nband * 2;
    
	gauss_n_max = 512;		        // default value: 2048
	gauss_n_min = 64;		        // default value: 64
    iter_max = 999;                 // default value: 999

    dlt_freq = 0.005;               
    eta_freq = 0.01;             
    freq_upp = 5.;                      
    freq_low = -5.;

    unit_omg = 0.01;
}

void Prmtr::set_values() {
    //model related
    hubbU = 8.;
    mu = 0.;
    // t.reset(4, 0.);
    // t[0] = 0. - mu; // onset energy.
    // t[1] = 1.;
    // t[2] =-0.3;
    // t[3] = 0.2;

    // fitting related
    fit_pow = 2.; // default value: 2.
    fit_rsd = 2.; // default value: 2.

    // NORG parameter.
    templet_restrain = {0, -0, -4,  0,  4,  0};
    templet_control =  {8, 20,  8,  8,  8, 20};
    ndiv = templet_control.size();
    norg_sets = 1;                                      // default value: 1
    nI2B = SUM(templet_control);                        // default value:
    iter_max_norg = 99;                                 // default
    // nooc_mode = STR("nooc");
    // nooc_mode = STR("cpnooc");
    nooc_mode = STR("cnooc");
    after_modify_prmtr();
    stage2_restrain = {0, -0, -3,  0,  3,  0};
    stage2_control =  {8, 20,  8,  8,  8, 20};
}

// we set first divison as impurity. The maximum number of cavity("-"); mean electron("+").
void Prmtr::after_modify_prmtr() const
{
    nI2B.reset(norg_sets, 0.);
    control_divs.reset(norg_sets + 1, ndiv);
    control_divs[0] = templet_restrain;
    for_Int(i, 0, norg_sets) control_divs[i + 1] = templet_control;
    for_Int(i, 0, norg_sets) nI2B[i] = SUM(control_divs[i + 1]) - control_divs[i + 1][0];
    nbath = SUM(nI2B);
    norbit = SUM(control_divs.tr()[0]) + nbath;
    npartical.reset(norg_sets, 0);
    for_Int(i, 0, norg_sets) npartical[i] = SUM(control_divs[i + 1])/2 ;
}

void Prmtr::change_the_norg_restrain_and_div(VecInt new_restrain, VecInt new_control) const {
    templet_restrain = new_restrain;
    templet_control = new_control;
    // nooc_mode = STR("nooc");
    ndiv = new_control.size();
    after_modify_prmtr();
}

void Prmtr::derive() {

    nfreq = 1 + Int_ROUND((freq_upp - freq_low) / dlt_freq);
    Re_z.reset(nfreq);
    for_Int(n, 0, nfreq) Re_z[n] = Cmplx(freq_low + n * dlt_freq, eta_freq);

    // max_omg = 4 * SQRT(SQR(hubbU) + DOT(t, t));    
    // max_omg = 2 * (ABS(hubbU) + 8 * SQRT(DOT(t, t) - t[0] * t[0]));
    max_omg = 2 * bandw;

    num_omg = Int_ROUND(max_omg / unit_omg / 2);
    Im_z.reset(num_omg);
	for_Int(n, 0, num_omg) Im_z[n] = Cmplx(0., (2 * n + 1) * unit_omg);


    fit_max_omg = max_omg/ 4. ;
    fit_num_omg = Int_ROUND(fit_max_omg / unit_omg / 2); // The default value: Int_ROUND(fit_max_omg / unit_omg / 2) change for speed reason.

    Str key_prmtr_val = STR("prmtr");
    ofx = iox + key_prmtr_val;


}

void Prmtr::print(std::ostream &os) const {
#define prmtr_print(var, comment) print(os, NAME(var), STR(var), comment)

    using namespace std;
    MatInt stage2(concat(stage2_restrain,stage2_control).mat(2,stage2_control.size()));    
	Str cnooc = nooc_mode;

    os << "// prmtr print bgn  " << present() << endl;
    prmtr_print(np, "number of mpi processes");
    prmtr_print(nbath, "number of bath orbital");

    prmtr_print(norbs, "number of spin-orbitals");
   	prmtr_print(cnooc, "Correlation nature orbital occupation constraint.");
    prmtr_print(control_divs, "to set the number of division and the shortcut restratin.");
    prmtr_print(stage2, "to set the number of division and the shortcut restratin @ stage2.");
    prmtr_print(hubbU, "The hubbard U interaction strength");
    // for_Int(i,0,t.size()) prmtr_print(t[i], "quarter bandwidth");
    prmtr_print(max_omg, "max_omg of imgreen");
    prmtr_print(eta_freq, "eta of omega + i * eta");
    prmtr_print(dlt_freq, "The real x-axis unit for retarded green function");
    prmtr_print(fit_max_omg, "");
    prmtr_print(fit_num_omg, "");
    prmtr_print(unit_omg*2, "The image x-axis unit for Matsubara Green's function");
    prmtr_print(imp_backup, "if you have the impurity's back up?");
    prmtr_print(nooc_mode, "the tight mode is refer to the correlation nature orbital constraint");
    prmtr_print(ofx, "output filename prefix");

    os << "// prmtr print end  " << present() << endl;

#undef prmtr_print
}