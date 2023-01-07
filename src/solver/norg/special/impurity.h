#pragma once

/*
coded by Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2021-02-25
*/

#include "specs.h"
#include "prmtr.h"
#include "green.h"
#include "bath.h"
#include "model.h"

class Impurity
{
private:
public:
  const MyMpi &mm; // parameters
  const Prmtr &p;  // parameters
  const Bath &bth;
  const Model &mdl;

  const Int ni; // number of impurity
  const Int ns; // number of sites,ns=ni+nb
  const Int nb; // number of bath sites
  MatReal V;
  MatReal E;
  MatReal h0; // hopping factors
private:
  // hopping  factors;when bath parameters is unusual,we just need to modify bath.hop()
  void set_factor();
  // void find_ve();
  // void sort_v_and_e();
  
public:
  Impurity(const MyMpi &mm_i, const Prmtr &prmtr_i, const Bath &bth_i, const Model &mdl, const Str &file = empty_str);
  void find_g0(Green &g0) const;
  void find_hb(Green &hb) const;
  void find_ve_for_test();
  void update()
  {
    // find_ve();
    // set_factor();
    // if(mm) WRN(NAV(h0));
  }
  void save() const;
};
