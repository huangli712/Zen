/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2013 - 2017
*/

#ifndef _TOOLBOX_H_
#define _TOOLBOX_H_

#include "stdpfx.h"
#include "vec.h"
#include "bijection.h"
#include "mat.h"
#include "number.h"
#include "bnmlfast.h"
#include "gauss.h"
#include "gaussj.h"
#include "fitmrq.h"
#include "tensor.h"
#include "iofmt.h"
#include "percent.h"
#include "nund.h"
#include "comdivs.h"
#include "onb.h"
#include "lanczos.h"
#include "bisection.h"
#include "slctsort.h"
//#include "uurandom.h"
#include "random.h"
#include "llintdouble.h"
#include "mymkl.h"
#include "mympi.h"
//#include "prng.h"
#include "myomp.h"
#include "ompidx.h"
#include "stat.h"
#include "sparsemat.h"
#include "msolver.h"


// Pauli arrays to define Pauli matrices
extern Cmplx pauli_array_s0[4];
extern Cmplx pauli_array_sx[4];
extern Cmplx pauli_array_sy[4];
extern Cmplx pauli_array_sz[4];

// Pauli matrices
extern const MatCmplx pauli_0;
extern const MatCmplx pauli_x;
extern const MatCmplx pauli_y;
extern const MatCmplx pauli_z;

inline Real entanglement_entropy(const VecReal &es)
{
#ifdef _ASSERTION_
	if (es.size() < 1) ERR(STR("es.size() < 1"));
#endif
	const Real eps = 100 * sqrt(Real(es.size())) * eps_Real;
	VecReal sp(es.size());
	sp = es;
	Real sum = SUM(es);
	if (ABS(sum - 1.) > eps) WRN("es not normalized, " + NAV2(sum - 1., eps));
	sp *= (1. / sum);
	std::sort(sp.p(), sp.p() + sp.size());
	sum = 0.;
	for_Int (i, 0, sp.size()) {
		const Real &v = sp[i];
		sum += v > 0. ? -v * log(v) : 0.;
	}
	return sum;
}

inline Real entanglement_entropy(const Vec<VecReal> &es)
{
	Idx n = 0;
	for_Idx (i, 0, es.size()) n += es[i].size();
	VecReal sp(n);
	n = 0;
	for_Idx (i, 0, es.size()) {
		for_Idx (j, 0, es[i].size()) sp[n++] = es[i][j];
	}
	return entanglement_entropy(sp);
}

#endif /* _TOOLBOX_H_ */
