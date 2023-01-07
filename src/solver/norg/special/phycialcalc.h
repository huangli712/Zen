#pragma once
/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

#include "specs.h"
#include "prmtr.h"
#include "nocspace.h"
#include "state.h"

// The here we get the sate from the NORG process, and we do some process to analyse the target state.
class PhycialCalc
{
private:



public:
    const MyMpi &mm;      // parameters
    const Prmtr &p;       // parameters
    const NocSpace &scsp; // NocSpace

    // const VecReal &state; // the target state
private:
    VecReal find_n_i(VecReal& state);
    VecReal find_Sp_i(VecReal& state);
    VecReal find_Sm_i(VecReal& state);
public:

PhycialCalc(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& s_i);


};