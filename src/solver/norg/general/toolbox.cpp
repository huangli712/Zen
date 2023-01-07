/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2013 - 2017
*/

#include "toolbox.h"

// Pauli arrays to define Pauli matrices
extern Cmplx pauli_array_s0[4] = { 1.,  0.,  0.,  1.};
extern Cmplx pauli_array_sx[4] = { 0.,  1.,  1.,  0.};
extern Cmplx pauli_array_sy[4] = { 0., -I,   I ,  1.};
extern Cmplx pauli_array_sz[4] = { 1.,  0.,  0., -1.};

// Pauli matrices
// extern const MatCmplx pauli_0(2, 2, pauli_array_s0);
// extern const MatCmplx pauli_x(2, 2, pauli_array_sx);
// extern const MatCmplx pauli_y(2, 2, pauli_array_sy);
// extern const MatCmplx pauli_z(2, 2, pauli_array_sz);
