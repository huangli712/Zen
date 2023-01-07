/*
code developed by
	Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2020-10-14

Copy from Rubin H. Landau, Manuel J. Paez, 
"Computational Phyiscs, Problem Solving With Computers", 
WILEY-INTERSCIENCE, 1997.
by Rong-Qiang He
Department of physics, Zhejiang University
Feb 26, 2009

Revise to F90 version from C version on Dec 25, 2010.
Then from F90 version to CPP version on 2020-10-14.

!********************************************************************
!* Points and weights for Gaussian quadrature                       *
!*       n	number of points                                        *
!* job = 0	rescaling uniformly between (a,b)                       *
!*       1	for integral (0,b) with 50% points inside (0,ab/(a+b))  *
!*       2	for integral (a,inf) with 50% inside (a,b+2a)           *
!********************************************************************

If you want to integrate a function f(x) in [a, b], the result is
SUM_{i=1,n} f(x[i]) * w[i];
with job = 0.
*/

#ifndef _GAUSS_H_
#define _GAUSS_H_

#include "stdpfx.h"
#include "vec.h"

class Gauss {
public:
	Int n;
	Int job;
	Real a;
	Real b;
	VecReal x;
	VecReal w;
public:
	inline Gauss(Int n_i, Int job_i, Real a_i, Real b_i);
private:
	void gauss_generate_abscissas_and_weights();
};

Gauss::Gauss(Int n_i, Int job_i, Real a_i, Real b_i) : n(n_i), job(job_i), a(a_i), b(b_i), x(n), w(n) {
	gauss_generate_abscissas_and_weights();
}

Real circle_function_for_gauss_integral_test(Real radius, Real x);
Real find_circle_area_for_gauss_integral_test(Real radius);

#endif /* _GAUSS_H_ */
