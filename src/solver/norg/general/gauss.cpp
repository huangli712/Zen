/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2020-10-15
*/

#include "gauss.h"

void Gauss::gauss_generate_abscissas_and_weights()
{
	if (n < 1 || n > 20000000 || job < 0 || job > 2) {
		ERR(NAV2(n, job));
	}
	const Real pi = 3.1415926535897932384626433832795;
	const Real eps = 5.e-15;

	const Int m = (n + 1) / 2;
	for_Int(i, 0, m) {
		Real pp;
		Real t = std::cos(pi * (i + 0.75) / (n + 0.5));
		Real t1 = 1.;
		while ((ABS(t - t1)) >= eps) {
			Real p1 = 1.;
			Real p2 = 0.;
			for_Int(j, 0, n) {
				Real p3 = p2;
				p2 = p1;
				p1 = ((2 * j + 1) * t * p2 - j * p3) / (j + 1);
			}
			pp = n * (t * p1 - p2) / (t * t - 1);
			t1 = t;
			t = t1 - p1 / pp;
		}
		x[i] = -t;
		x[n - i - 1] = t;
		w[i] = 2. / ((1 - t * t) * pp * pp);
		w[n - i - 1] = w[i];
	}
	
	if (job == 0) {
		const Real half = (b - a) / 2.;
		const Real midpoint = (b + a) / 2.;
		for_Int(i, 0, n) {
			x[i] = x[i] * half + midpoint;
			w[i] = w[i] * half;
		}
	}
	
	if (job == 1) {
		for_Int(i, 0, n) {
			Real tmp = b + a - (b - a) * x[i];
			x[i] = a * b * (1 + x[i]) / tmp;
			w[i] = w[i] * 2 * a * b * b / (tmp * tmp);
		}
	}

	if (job == 2) {
		for_Int(i, 0, n) {
			Real tmp = 1 - x[i];
			x[i] = (b * x[i] + b + a + a) / tmp;
			w[i] = w[i] * 2 * (a + b) / (tmp * tmp);
		}
	}
}

Real circle_function_for_gauss_integral_test(Real radius, Real x) {
	return SQRT(radius * radius - x * x);
}

Real find_circle_area_for_gauss_integral_test(Real radius)
{
	Real f(Real x);

	const Real pi = 3.1415926535897932384626433832795;
	Int gauss_n = 1000;
	Int gauss_job = 0;
	Real a = 0.;
	Real b = radius;
	Gauss gauss(gauss_n, gauss_job, a, b);
	Real sum = 0.;
	for_Int(i, 0, gauss.n) {
		sum += circle_function_for_gauss_integral_test(radius, gauss.x[i]) * gauss.w[i];
	}
	return sum * 4;
}
