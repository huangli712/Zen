#pragma once


/*
coded by Yi - Heng Tian(ructyh@ruc.edu.cn, RUC, China)
date 20211012
*/

#include "specs.h"
#include "prmtr.h"
#include "model.h"
#include "green.h"

class Phy {
public:
	const MyMpi& mm;
	const Prmtr& p;
	const Model& mdl;
    StdVecCmplx vec_k;

  private:
	void set_k();
    MatCmplx find_dk(Real kx, Real ky);
    MatCmplx find_pk(Real kx, Real ky);
    //Cmplx find_pk(Real kx, Real ky);

  public:
    Phy(const MyMpi& mm_i, const Prmtr& prmtr_i, const Model& mdl_i);
	//find_spectral function
    Real find_M();
    Cmplx find_D();
    Cmplx find_P();
};