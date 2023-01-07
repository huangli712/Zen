#pragma once

/*
code	by	Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China) date 2013 - 2017
modify	by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

#include "specs.h"
#include "prmtr.h"
#include "green.h"

class HybErr {
public:
	const Prmtr& p;			// parameters
	const ImGreen& hb;		// hybridization function
	const Int nw;			// nw = p.fit_num_omg
	Int nb;			// nb = p.nbath
	VecInt x;
	VecReal y;
	VecReal sig;
	Int mode;
public:
	HybErr(const Prmtr& p_i, const ImGreen& hb_i);

	HybErr(const Prmtr& p_i, const ImGreen& hb_i, const Int& nb_i, Int orb_i,Int mode);

	// return y = f(x; a);
	Real operator()(const Int x, const VecReal& a) const {
		VecReal temp;
		if (x < 2 * nw) {
			const Int n = x < nw ? x : x - nw;
			const Real omgn = p.Imomg(n);
			if (mode == 1)
			{
				Int nb_tmp = nb / 2;
				if (nb % 2 == 0)
				{
					// const VecCmplx E = cmplx(concat(VecReal(nb_tmp, a.p()), -VecReal(nb_tmp, a.p())));
					temp.sm(nb_tmp, a.p()); const VecCmplx E = cmplx(concat(temp, -temp));
					const VecCmplx S = INV(E - VecCmplx(2 * nb_tmp, I * omgn));
					temp.sm(nb_tmp, a.p() + nb_tmp); const VecCmplx V = cmplx(concat(temp, temp));
					const VecCmplx Vco = V.co();
					Cmplx hyb = DOT(Vco, S * Vco);
					return x < nw ? real(hyb) : imag(hyb);
				}
				else
				{
					if (nb == 1)
					{
						const VecCmplx V = cmplx(a);
						const VecCmplx S = -INV(VecCmplx(nb, I * omgn));
						const VecCmplx Vco = V.co();
						Cmplx hyb = DOT(Vco, S * Vco);
						return x < nw ? real(hyb) : imag(hyb);
					}
					else if (nb > 1)
					{
						VecReal E_tmp; E_tmp.sm(nb_tmp, a.p());
						VecReal V_tmp_left;  V_tmp_left.sm(nb_tmp + 1, a.p() + nb_tmp);
						VecReal V_tmp_right;  V_tmp_right.sm(nb_tmp, a.p() + nb_tmp);

						VecReal tmp0(1, 0.);
						const VecCmplx E = cmplx(concat(concat(E_tmp, tmp0), -E_tmp));
						const VecCmplx S = INV(E - VecCmplx(2 * nb_tmp + 1, I * omgn));
						const VecCmplx V = cmplx(concat(V_tmp_left, V_tmp_right));
						const VecCmplx Vco = V.co();
						Cmplx hyb = DOT(Vco, S * Vco);
						return x < nw ? real(hyb) : imag(hyb);
					}
				}
			}
			else if (mode == 0)
			{
				const VecCmplx E = cmplx(temp.sm(nb, a.p()));
				const VecCmplx S = INV(E - VecCmplx(nb, I * omgn));
				const VecCmplx V = cmplx(temp.sm(nb, a.p() + nb));
				const VecCmplx Vco = V.co();
				Cmplx hyb = DOT(Vco, S * Vco);
				return x < nw ? real(hyb) : imag(hyb);
			}
		}
		else if (x == 2 * nw + 0) {

			if (mode == 1)
			{
				Int nb_tmp = nb / 2;
				if (nb % 2 == 0)
				{
					const VecReal E = concat(temp.sm(nb_tmp, a.p()), -temp.sm(nb_tmp, a.p()));
					const VecReal E2 = E * E;
					return DOT(E2, E2);
				}
				else
				{
					if (nb == 1)
					{
						const VecReal V = a;
						const VecReal Vco = V.co();
						Real hyb = DOT(Vco, Vco);

						
						return hyb;
						
					}
					else if (nb > 1)
					{
						const VecReal E_tmp = temp.sm(nb_tmp, a.p());
						VecReal tmp0(1, 0.);
						const VecReal E = concat(concat(E_tmp, tmp0), -E_tmp);
						const VecReal E2 = E * E;
						return DOT(E2, E2);
					}
					
				}
			}
			else if (mode == 0)
			{
				const VecReal E = temp.sm(nb, a.p());
				const VecReal E2 = E * E;
				return DOT(E2, E2);
			}

			
		}
		else if (x == 2 * nw + 1) {
			if (mode == 1)
			{
				Int nb_tmp = nb / 2;
				if (nb % 2 == 0)
				{
					const VecReal V = concat(temp.sm(nb_tmp, a.p() + nb_tmp), temp.sm(nb_tmp, a.p() + nb_tmp));
					const VecReal Vco = V.co();
					Real hyb = DOT(Vco, Vco);
					return hyb;
				}
				else
				{
					if (nb == 1)
					{
						const VecReal V = a;
						const VecReal Vco = V.co();
						Real hyb = DOT(Vco, Vco);

						
						return hyb;
					}
					else if (nb > 1)
					{
						const VecReal V_tmp_left = temp.sm(nb_tmp + 1, a.p() + nb_tmp);
						const VecReal V_tmp_right = temp.sm(nb_tmp, a.p() + nb_tmp);
						const VecReal V = concat(V_tmp_left, V_tmp_right);
						const VecReal Vco = V.co();
						Real hyb = DOT(Vco, Vco);
						return hyb;
					}
					
				}
			}
			else if (mode == 0)
			{
				const VecReal V = temp.sm(nb, a.p() + nb);
				const VecReal Vco = V.co();
				Real hyb = DOT(Vco, Vco);
				return hyb;
			}
			
		}
	}
	// return y = f(x; a) and dy/da
	void operator()(const Int x, const VecReal& a, Real& y, VecReal& dyda) const;
	// return relative_err = SQRT(chi_sqr)
	Real operator()(const VecReal& a) const {
		VecReal fx(x.size());
		for_Int(i, 0, x.size()) {
			fx[i] = (*this)(x[i], a);
		}
		VecReal relative_dev = (fx - y) * INV(sig);
		return SQRT(DOT(relative_dev, relative_dev));
	}
	// return relative_err = SQRT(the curve part of chi_sqr)
	Real err_curve(const VecReal& a) const {
		VecReal fx(x.size());
		for_Int(i, 0, x.size()) {
			fx[i] = (*this)(x[i], a);
		}
		//fx[2 * nw] = y[2 * nw];
		//fx[2 * nw + 1] = y[2 * nw + 1];
		VecReal relative_dev = (fx - y) * INV(sig);
		return SQRT(DOT(relative_dev, relative_dev));
	}
	// return relative_err = SQRT(the part of ose regularization)
	Real err_reg(const VecReal& a) const {
		const Int i = 2 * nw + 0;
		Real fx = (*this)(x[i], a);
		Real relative_dev = (fx - y[i]) * INV(sig[i]);
		return SQRT(relative_dev * relative_dev);
	}
	// return relative_err = SQRT(the part of bath sum rule)
	Real err_bsr(const VecReal& a) const {
		const Int i = 2 * nw + 1;
		Real fx = (*this)(x[i], a);
		Real relative_dev = (fx - y[i]) * INV(sig[i]);
		return SQRT(relative_dev * relative_dev);
	}
	void write_xysig(Int iter_cnt) const {
		OFS ofs(iox + "zic" + prefill0(iter_cnt, 3) + ".xysig.txt");
		using namespace std;
		ofs << setw(6) << "x" << "  " << setw(w_Real) << "y" << "  " << setw(w_Real) << "sig" << endl;
		ofs << iofmt("sci");
		for_Int(i, 0, x.size()) {
			ofs << setw(6) << x[i] << "  " << setw(w_Real) << y[i] << "  " << setw(w_Real) << sig[i] << endl;
		}
	}
};