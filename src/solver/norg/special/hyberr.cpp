/*
code	by	Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China) date 2013 - 2017
modify	by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/
#include "hyberr.h"

HybErr::HybErr(const Prmtr& p_i, const ImGreen& hb_i) :
	p(p_i), hb(hb_i), nw(p.fit_num_omg), nb(p.nbath),
	x(2 * nw + 1), y(2 * nw + 1), sig(2 * nw + 1)
{
	// curve
	{
		Real mag_real = 0.;
		Real mag_imag = 0.;
		VecReal wght(2 * nw);
		for_Int(n, 0, nw) {
			x[n] = n;
			y[n] = real(hb[n][0][0]);
			mag_real += SQR(y[n]);
			x[nw + n] = nw + n;
			y[nw + n] = imag(hb[n][0][0]);
			mag_imag += SQR(y[nw + n]);
			// wght[n] = wght[nw + n] = std::pow(p.mbomg(n) + p.fit_rsd, -p.fit_pow);
			wght[n] = wght[nw + n] = 1;
		}
		mag_real = SQRT(1 + mag_real / nw);
		mag_imag = SQRT(1 + mag_imag / nw);
		wght *= INV(SUM(wght));
		for_Int(n, 0, nw) {
			sig[n] = mag_real / SQRT(wght[n]);
			sig[nw + n] = mag_imag / SQRT(wght[nw + n]);
		}
	}
	// the part of ose regularization
	if(x.size()>= 2 * nw + 1){
		x[2 * nw + 0]	= 2 * nw + 0;
		y[2 * nw + 0]	= 0.;
		sig[2 * nw + 0] = nb * std::pow(256 * p.fit_max_omg, 4);
	}
	// the part of bath sum rule
	// if (x.size() >= 2 * nw + 2) {
	// 	x[2 * nw + 1]	= 2 * nw + 1;
	// 	y[2 * nw + 1]	= 0.;
	// 	sig[2 * nw] = 0.1 * (1 + ABS(p.bsr));  //Change to the part of ose regularization
	// }
}

HybErr::HybErr(const Prmtr& p_i, const ImGreen& hb_i, const Int& nb_i, Int orb_i, Int mode_i) :
	p(p_i), hb(hb_i), nw(p.fit_num_omg), nb(nb_i), mode(mode_i),
	x(2 * nw + 2), y(2 * nw + 2), sig(2 * nw + 2)
{
	// curve
	{
		
		Real mag_real = 0.;
		Real mag_imag = 0.;
		VecReal wght(2 * nw);

		for_Int(n, 0, nw) {
			x[n] = n;
			y[n] = real(hb[n][0][0]);
			mag_real += SQR(y[n]);
			x[nw + n] = nw + n;
			y[nw + n] = imag(hb[n][0][0]);
			mag_imag += SQR(y[nw + n]);
			// wght[n] = wght[nw + n] = std::pow(p.mbomg(n) + p.fit_rsd, -p.fit_pow);
			wght[n] = wght[nw + n] = 1;
		}

		mag_real = SQRT(1 + mag_real / nw);
		mag_imag = SQRT(1 + mag_imag / nw);
		wght *= INV(SUM(wght));
		for_Int(n, 0, nw) {
			sig[n] = mag_real / SQRT(wght[n]);
			sig[nw + n] = mag_imag / SQRT(wght[nw + n]);
		}
	}
	

	
	// the part of ose regularization
	x[2 * nw + 0] = 2 * nw + 0;
	y[2 * nw + 0] = 0.;
	if (mode == 0)
	{
		sig[2 * nw + 0] = nb * std::pow(256 * p.fit_max_omg, 4);
	}
	/* 
	else if (mode == 1)
	{
		if (nb == 1)
		{
			//sig[2 * nw + 0] = nb * std::pow(8 * p.fit_max_omg, 4);
			x[2 * nw + 0] = 2 * nw + 0;
			y[2 * nw + 0] = p.t[orb_i] * p.t[orb_i];
    			sig[2 * nw + 0] = 0.1 * (1 + p.t[orb_i] * p.t[orb_i]);
		}
		else if (nb > 1)
		{
			sig[2 * nw + 0] = (nb / 2) * std::pow(256 * p.fit_max_omg, 4);
		}
		
	}
	 */
	
	/* 
	// the part of bath sum rule
	x[2 * nw + 1] = 2 * nw + 1;
	y[2 * nw + 1] = p.t[orb_i] * p.t[orb_i];
    // sig[2 * nw + 1] = 0.1 * (1 + p.t[orb_i] * p.t[orb_i]);
    sig[2 * nw + 1] = 50 * 0.1 * (1 + p.t[orb_i] * p.t[orb_i]);
     */
}

/*
chi_sqr = sum_i [(y_i - hb(x_i)) / sig_i]^2 + [(VV^+ - bsr) / sig_bsr]^2 + ose_regularization
chi_sqr = sum_i weight_i * [(y_i - hb(x_i)) / mag] ^ 2
        + weight_bsr * [(VV^+ - bsr) / mag_bsr] ^ 2
		+ ose_regularization
E = ose;			// on-site energies for bath sites
V = hop;			// hopping coefficients between impurity site and bath sites
S = INV(diag_matrix(E - i * omgn));
hyb = hb(omgn) = V * S * V^+ = SUM_i V_i * S_{i,i} * V_i^*
D = d(hb) / d(a) = concat(d(hb) / d(E), d(hb) / d(V))
d(S_{i,i}) / d(E_i) = d(1/(E_i - i * omgn)) / d(E_i) = -(E_i - i * omgn)^{-2} = -S{i,i}^2
(D_E)_i = d(hyb)/d(E_i) = V_i * d(S_{i,i}) / d(E_i) * V_i^* = V_i * -S{i,i}^2 * V_i^*
*/
void HybErr::operator()(const Int x, const VecReal& a, Real& y, VecReal& dyda) const
{
	VecReal temp;
	if (x < 2 * nw) {
		const Int n = x < nw ? x : x - nw;
		const Cmplx iomgn = p.Imomg(n);
		if (mode == 1)
		{
			Int nb_tmp = nb / 2;
			if (nb % 2 == 0)
			{
				const VecReal E_tmp = temp.sm(nb_tmp, a.p());
				const VecReal V_tmp = temp.sm(nb_tmp, a.p() + nb_tmp);

				const VecCmplx E = cmplx(concat(E_tmp, -E_tmp));
				const VecCmplx S = INV(E - VecCmplx(2 * nb_tmp, iomgn));
				const VecCmplx V = cmplx(concat(V_tmp, V_tmp));
				const VecCmplx Vco = V.co();
				Cmplx hyb = DOT(Vco, S * Vco);
				y = x < nw ? real(hyb) : imag(hyb);
				//VecCmplx D_E = V * Cmplx(-1.) * S * S * Vco;

				VecCmplx D_E = cmplx(V_tmp) * INV(cmplx(E_tmp) + VecCmplx(nb_tmp, iomgn)) * INV(cmplx(E_tmp) + VecCmplx(nb_tmp, iomgn)) * (cmplx(V_tmp.co())) - cmplx(V_tmp) * INV(cmplx(E_tmp) - VecCmplx(nb_tmp, iomgn)) * INV(cmplx(E_tmp) - VecCmplx(nb_tmp, iomgn)) * (cmplx(V_tmp.co()));
				//VecCmplx D_V = V * S + S * Vco;
				VecCmplx D_V = cmplx(V_tmp) * (INV(cmplx(E_tmp) - VecCmplx(nb_tmp, iomgn)) - INV(cmplx(E_tmp) + VecCmplx(nb_tmp, iomgn))) + (INV(cmplx(E_tmp) - VecCmplx(nb_tmp, iomgn)) - INV(cmplx(E_tmp) + VecCmplx(nb_tmp, iomgn))) * (cmplx(V_tmp.co()));
				VecCmplx D = concat(D_E, D_V);
				dyda = x < nw ? real(D) : imag(D);
			}
			else
			{
				if (nb == 1)
				{
					const VecCmplx V = cmplx(a);
					const VecCmplx S = -INV(VecCmplx(nb, iomgn));
					const VecCmplx Vco = V.co();
					Cmplx hyb = DOT(Vco, S * Vco);
					y = x < nw ? real(hyb) : imag(hyb);
					VecCmplx D_V = V * (INV(- VecCmplx(nb, iomgn))) + (INV(- VecCmplx(nb, iomgn))) * (V.co());
					dyda = x < nw ? real(D_V) : imag(D_V);
				}
				else if (nb > 1)
				{
					const VecReal E_tmp = temp.sm(nb_tmp, a.p());
					const VecReal V_tmp_left = temp.sm(nb_tmp + 1, a.p() + nb_tmp);
					const VecReal V_tmp_right = temp.sm(nb_tmp, a.p() + nb_tmp);

					VecReal tmp0(1, 0.);
					const VecCmplx E = cmplx(concat(concat(E_tmp, tmp0), -E_tmp));
					const VecCmplx S = INV(E - VecCmplx(2 * nb_tmp + 1, iomgn));
					const VecCmplx V = cmplx(concat(V_tmp_left, V_tmp_right));
					const VecCmplx Vco = V.co();
					Cmplx hyb = DOT(Vco, S * Vco);
					y = x < nw ? real(hyb) : imag(hyb);
					//VecCmplx D_E = V * Cmplx(-1.) * S * S * Vco;

					VecCmplx D_E = cmplx(V_tmp_right) * INV(cmplx(E_tmp) + VecCmplx(nb_tmp, iomgn)) * INV(cmplx(E_tmp) + VecCmplx(nb_tmp, iomgn)) * (cmplx(V_tmp_right.co())) - cmplx(V_tmp_right) * INV(cmplx(E_tmp) - VecCmplx(nb_tmp, iomgn)) * INV(cmplx(E_tmp) - VecCmplx(nb_tmp, iomgn)) * (cmplx(V_tmp_right.co()));
					//VecCmplx D_V = V * S + S * Vco;
					VecCmplx D_V_tmp(1, -V_tmp_left[nb_tmp] * cmplx(1.) / iomgn - (cmplx(1.) / iomgn) * cnjg(V_tmp_left[nb_tmp]));
					VecCmplx D_V = concat(cmplx(V_tmp_right) * (INV(cmplx(E_tmp) - VecCmplx(nb_tmp, iomgn)) - INV(cmplx(E_tmp) + VecCmplx(nb_tmp, iomgn))) + (INV(cmplx(E_tmp) - VecCmplx(nb_tmp, iomgn)) - INV(cmplx(E_tmp) + VecCmplx(nb_tmp, iomgn))) * (cmplx(V_tmp_right.co())), D_V_tmp);
					VecCmplx D = concat(D_E, D_V);
					dyda = x < nw ? real(D) : imag(D);
				}
				
			}
		}
		else if (mode == 0)
		{
			const VecCmplx E = cmplx(temp.sm(nb, a.p()));
			const VecCmplx S = INV(E - VecCmplx(nb, iomgn));
			const VecCmplx V = cmplx(temp.sm(nb, a.p() + nb));
			const VecCmplx Vco = V.co();
			Cmplx hyb = DOT(Vco, S * Vco);
			y = x < nw ? real(hyb) : imag(hyb);
			VecCmplx D_E = V * Cmplx(-1.) * S * S * Vco;
			VecCmplx D_V = V * S + S * Vco;
			VecCmplx D = concat(D_E, D_V);
			dyda = x < nw ? real(D) : imag(D);
		}
	}
	else if (x == 2 * nw + 0) { // the part of ose regularization
		if (mode == 1)
		{
			Int nb_tmp = nb / 2;
			if (nb % 2 == 0)
			{
				const VecReal E_tmp = temp.sm(nb_tmp, a.p());
				const VecReal E = concat(E_tmp, -E_tmp);
				const VecReal E2 = E * E;
				y = DOT(E2, E2);
				VecReal D_E = 8. * (E_tmp * E_tmp * E_tmp);
				VecReal D_V = VecReal(nb_tmp, 0.);
				VecReal D = concat(D_E, D_V);
				dyda = D;
			}
			else
			{
				if (nb == 1)
				{
					const VecReal V = a;
					const VecReal Vco = V.co();
					Real hyb = DOT(Vco, Vco);
					y = hyb;
					VecReal D_V = 2. * V;
					dyda = D_V;
				}
				else if (nb > 1)
				{
					const VecReal E_tmp = temp.sm(nb_tmp, a.p());
					VecReal tmp0(1, 0.);
					const VecReal E = concat(concat(E_tmp, tmp0), -E_tmp);
					const VecReal E2 = E * E;
					y = DOT(E2, E2);
					VecReal D_E = 8. * (E_tmp * E_tmp * E_tmp);
					VecReal D_V = VecReal(nb_tmp + 1, 0.);
					VecReal D = concat(D_E, D_V);
					dyda = D;
				}
			}
		}
		else if (mode == 0)
		{
			const VecReal E = temp.sm(nb, a.p());
			const VecReal E2 = E * E;
			y = DOT(E2, E2);
			VecReal D_E = 4. * E2 * E;
			VecReal D_V = VecReal(nb, 0.);
			VecReal D = concat(D_E, D_V);
			dyda = D;
		}
	}
	else if (x == 2 * nw + 1) { // the part of bath sum rule
		if (mode == 1)
		{
			Int nb_tmp = nb / 2;
			if (nb % 2 == 0)
			{
				const VecReal V_tmp = temp.sm(nb_tmp, a.p() + nb_tmp);
				const VecReal V = concat(V_tmp, V_tmp);
				const VecReal Vco = V.co();
				Real hyb = DOT(Vco, Vco);
				y = hyb;
				VecReal D_E = VecReal(nb_tmp, 0.);
				VecReal D_V = 2. * V_tmp + 2. * (V_tmp.co());
				VecReal D = concat(D_E, D_V);
				dyda = D;
			}
			else
			{
				if (nb == 1)
				{
					const VecReal V = a;
					const VecReal Vco = V.co();
					Real hyb = DOT(Vco, Vco);
					y = hyb;
					VecReal D_V = 2. * V;
					dyda = D_V;
				}
				else if (nb > 1)
				{
					const VecReal V_tmp_left = temp.sm(nb_tmp + 1, a.p() + nb_tmp);
					const VecReal V_tmp_right = temp.sm(nb_tmp, a.p() + nb_tmp);
					const VecReal V = concat(V_tmp_left, V_tmp_right);
					const VecReal Vco = V.co();
					Real hyb = DOT(Vco, Vco);
					y = hyb;
					VecReal D_E = VecReal(nb_tmp, 0.);
					VecReal D_V_tmp(1, V_tmp_left[nb_tmp] + cnjg(V_tmp_left[nb_tmp]));
					VecReal D_V = concat(2. * V_tmp_right + 2. * (V_tmp_right.co()), D_V_tmp);
					VecReal D = concat(D_E, D_V);
					dyda = D;
				}
			}
		}
		else if (mode == 0)
		{
			const VecReal V = temp.sm(nb, a.p() + nb);
			const VecReal Vco = V.co();
			Real hyb = DOT(Vco, Vco);
			y = hyb;
			VecReal D_E = VecReal(nb, 0.);
			VecReal D_V = V + Vco;
			VecReal D = concat(D_E, D_V);
			dyda = D;
		}

	}
	else {
		ERR(NAV(x));
	}
}