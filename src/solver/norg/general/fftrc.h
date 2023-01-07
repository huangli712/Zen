#pragma once

/*******************************************************************************
* Forward FFT from real to complex (with the conjugate-even symmetry)
* Z(k1, ..., kd) = SUM_{n1 = 0, ..., N1-1; ...; nd = 0, ..., Nd-1} * R(n1, ..., nd) *
* EXP(-2 * pi * i / N1 * k1 * n1) * ... * EXP(-2 * pi * i / Nd * kd * nd)
* DFTI_CONJUGATE_EVEN_STORAGE: storage scheme for a conjugate-even domain.
* DFTI_COMPLEX_REAL is deprecated; DFTI_COMPLEX_COMPLEX is used.
* Because the input sequence R is real-valued, the mathematical result Z has
* conjugate-even symmetry: Z(k1, k2, ..., kd) = conjugate (Z(N1-k1, N2-k2, ..., Nd-kd)).
* With the conjugate-even symmetry, approximately a half of the result suffices to fully reconstruct it.
* In the Intel(R) oneAPI Math Kernel Library FFT interface, the halved dimension is the last dimension.
* It suffices to store elements Z(k1, ..., kd) for the following indices:
* ki = 0, бн, Ni-1, where i = 1, ..., d-1,
* kd = 0, ..., [Nd/2].
* NOTE: only d = 2 is implemented below.
* Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China) 2021-03-12
*******************************************************************************/

#include "stdpfx.h"
#include "vec.h"
#include "mat.h"
#include "mkl_service.h"
#include "mkl_dfti.h"

// forward FFT from real to complex(with the conjugate - even symmetry)
// Z(k1, ..., kd) = SUM_{ n1 = 0, ..., N1 - 1; ...; nd = 0, ..., Nd - 1 } * R(n1, ..., nd)
// * EXP(-2 * pi * i / N1 * k1 * n1) * ... * EXP(-2 * pi * i / Nd * kd * nd)
// x is a matrix, which means d = 2
MatCmplx fftrc(const MatReal& x);

// forward FFT from real to complex(with the conjugate - even symmetry)
// Z(k1, ..., kd) = SUM_{ n1 = 0, ..., N1 - 1; ...; nd = 0, ..., Nd - 1 } * R(n1, ..., nd)
// * EXP(-2 * pi * i / N1 * k1 * n1) * ... * EXP(-2 * pi * i / Nd * kd * nd)
// x is a vector, which means d = 1
VecCmplx fftrc(const VecReal& x);
