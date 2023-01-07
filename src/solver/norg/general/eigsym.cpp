/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2012 - 2017
*/

#include "eigsym.h"
#include "slctsort.h"

// the eigenpairs are in ascending order
// columns of z stores eigenvectors on output
// created by Rong-Qiang He, 2014-01-15
// ifv indicates if eigenvectors are requested, ifv has default value 1
// the return value is ierr, ierr = -1 indicates success
// ierr >= 0 indicates the number of eigenpairs have already been determined
// ierr-th eigenpair can not be determined in itr_max iterations
Int eigsym(MatReal &z, VecReal &d, Int ifv)
{
#ifdef _CHECK_DIMENSION_MATCH_
	if (z.nrows() != z.ncols() || z.nrows() != d.size())
		ERR("z.nrows() != z.ncols() || z.nrows() != d.size()"
		", z.nrows() = " + STR(z.nrows()) + 
		", z.ncols() = " + STR(z.ncols()) + 
		", d.size() = " + STR(d.size()));
#endif
	VecReal e(d.size());
	tred2(z, d, e, ifv);
	return tqli(z, d, e, ifv);
}

// transform a symmetric matrix z to tridiagonal
// if ifv == 1, z contains the transformation matrix on output
// e[0] = 0 on output
void tred2(MatReal &z, VecReal &d, VecReal &e, Int ifv)
{
#ifdef _CHECK_DIMENSION_MATCH_
	if (z.nrows() != z.ncols() || z.nrows() != d.size() || z.nrows() != e.size())
		ERR("z.nrows() != z.ncols() || z.nrows() != d.size() || z.nrows() != e.size()"
		", z.nrows() = " + STR(z.nrows()) + 
		", z.ncols() = " + STR(z.ncols()) + 
		", d.size() = " + STR(d.size()) + 
		", e.size() = " + STR(e.size()));
#endif
	Int n = (Int)d.size();
	Int l, k, j, i;
	Real scale, hh, h, g, f;
	if (n < 1) return;
	for (i = n - 1; i > 0; --i) {
		l = i - 1;
		h = scale = 0.;
		if (l > 0) {
			for (k = 0; k < i; k++) scale += ABS(z[i][k]);
			if (scale == 0.)
				e[i] = z[i][l];
			else {
				for (k = 0; k < i; k++) {
					z[i][k] /= scale;
					h += z[i][k] * z[i][k];
				}
				f = z[i][l];
				g = (f >= 0. ? -sqrt(h) : sqrt(h));
				e[i] = scale * g;
				h -= f * g;
				z[i][l] = f - g;
				f = 0.;
				for (j = 0; j < i; ++j) {
					if (ifv)
						z[j][i] = z[i][j] / h;
					g = 0.;
					for (k = 0; k < j + 1; k++) g += z[j][k] * z[i][k];
					for (k = j + 1; k < i; k++) g += z[k][j] * z[i][k];
					e[j] = g / h;
					f += e[j] * z[i][j];
				}
				hh = f / (h + h);
				for (j = 0; j < i; ++j) {
					f = z[i][j];
					e[j] = g = e[j] - hh * f;
					for (k = 0; k < j + 1; k++) z[j][k] -= (f * e[k] + g * z[i][k]);
				}
			}
		} else
			e[i] = z[i][l];
		d[i] = h;
	}
	if (ifv) d[0] = 0.;
	e[0] = 0.;
	for (i = 0; i < n; ++i) {
		if (ifv) {
			if (d[i] != 0.) {
				for (j = 0; j < i; ++j) {
					g = 0.;
					for (k = 0; k < i; k++) g += z[i][k] * z[k][j];
					for (k = 0; k < i; k++) z[k][j] -= g * z[k][i];
				}
			}
			d[i] = z[i][i];
			z[i][i] = 1.;
			for (j = 0; j < i; ++j) z[j][i] = z[i][j] = 0.;
		} else {
			d[i] = z[i][i];
		}
	}
}

// revised by Rong-Qiang He, 2014-01-15 from Numerical Recipies, Third Edition
// the eigenpairs are in ascending order
// columns of z stores eigenvectors on output
// e[0] = 0 on input
// only revision is commented here
// ifv indicates if eigenvectors are requested, ifv has default value 1
// if ifv == 1 and tqli is     preceded by a tred2, z on input should be the z of tred2 on output
// if ifv == 1 and tqli is not preceded by a tred2, z should be identity matrix on input
// the return value is ierr, ierr = -1 indicates success
// ierr  >=  0 indicates the number of eigenpairs have already been determined
// ierr-th eigenpair can not be determined in itr_max iterations
Int tqli(MatReal &z, VecReal &d, VecReal &e, Int ifv)
{
#ifdef _CHECK_DIMENSION_MATCH_
	if (z.nrows() != z.ncols() || z.nrows() != d.size() || z.nrows() != e.size())
		ERR("z.nrows() != z.ncols() || z.nrows() != d.size() || z.nrows() != e.size()"
		", z.nrows() = " + STR(z.nrows()) + 
		", z.ncols() = " + STR(z.ncols()) + 
		", d.size() = " + STR(d.size()) + 
		", e.size() = " + STR(e.size()));
#endif
	Int n = (Int)d.size();
	Int m, l, j, i, k, ierr;
	Int itr_max = Int(60 * log10(eps_Real) / log10(eps_double));
	Real s, r, p, g, f, dd, c, b;

	ierr = -1;
	if (n < 1) return ierr;
	for (i = 1; i < n; ++i) e[i - 1] = e[i];
	e[n - 1] = 0.;
	for (l = 0; ierr == -1 && l < n; ++l) {
		j = 0;
		do {
			for (m = l; m < n - 1; ++m) {
				dd = ABS(d[m]) + ABS(d[m + 1]);
				if (dd + ABS(e[m]) == dd) break;
			}
			if (m != l) {
				if (j++ == itr_max) {
					ierr = l;
					ERR("too many iterations in tqli"
						", itr_max = " + STR(itr_max) + 
						", j = " + STR(j) + 
						", m = " + STR(m) + 
						", l = " + STR(l) + 
						", z.nrows() = " + STR(z.nrows()));
					break;
				}
				g = (d[l + 1] - d[l]) / (2. * e[l]);
				r = pythag(g, 1.);
				g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
				s = c = 1.;
				p = 0.;
				for (i = m - 1; i >= l; --i) {
					f = s * e[i];
					b = c * e[i];
					e[i + 1] = (r = pythag(f, g));
					if (r == 0.) {
						d[i + 1] -= p;
						e[m] = 0.;
						break;
					}
					s = f / r;
					c = g / r;
					g = d[i + 1] - p;
					r = (d[i] - g) * s + 2. * c * b;
					d[i + 1] = g + (p = s * r);
					g = c * r - b;
					if (ifv) {
						for (k = 0; k < n; k++) {
							f = z[k][i + 1];
							z[k][i + 1] = s * z[k][i] + c * f;
							z[k][i] = c * z[k][i] - s * f;
						}
					}
				}
				if (r == 0. && i >= l) continue;
				d[l] -= p;
				e[l] = g;
				e[m] = 0.;
			}
		} while (m != l);
	}
	if (ifv)
		slctsort_matcol(d, z);
	else
		sort(d);
	return ierr;
}

// special version of tqli, z in tqlis coincides exactly with the last row of z in tqli
// revised by Rong-Qiang He, 2014-01-15 from Numerical Recipies, Third Edition
Int tqlis(VecReal &z, VecReal &d, VecReal &e, Int ifv)
{
#ifdef _CHECK_DIMENSION_MATCH_
	if (z.size() != d.size() || z.size() != e.size()) ERR("z.size() != d.size() || z.size() != e.size()"
		", z.size() = " + STR(z.size()) + 
		", d.size() = " + STR(d.size()) + 
		", e.size() = " + STR(e.size()));
#endif
	Int n = (Int)d.size();
	Int m, l, j, i, ierr;
	Int itr_max = Int(60 * log10(eps_Real) / log10(eps_double));
	Real s, r, p, g, f, dd, c, b;

	ierr = -1;
	if (n < 1) return ierr;
	for (i = 1; i < n; ++i) e[i - 1] = e[i];
	e[n - 1] = 0.;
	for (l = 0; ierr == -1 && l < n; ++l) {
		j = 0;
		do {
			for (m = l; m < n - 1; ++m) {
				dd = ABS(d[m]) + ABS(d[m + 1]);
				if (dd + ABS(e[m]) == dd) break;
			}
			if (m != l) {
				if (j++ == itr_max) {
					ierr = l;
					ERR("too many iterations in tqlis"
						", itr_max = " + STR(itr_max) + 
						", j = " + STR(j) + 
						", m = " + STR(m) + 
						", l = " + STR(l) + 
						", z.size() = " + STR(z.size()));
					break;
				}
				g = (d[l + 1] - d[l]) / (2. * e[l]);
				r = pythag(g, 1.);
				g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
				s = c = 1.;
				p = 0.;
				for (i = m - 1; i >= l; --i) {
					f = s * e[i];
					b = c * e[i];
					e[i + 1] = (r = pythag(f, g));
					if (r == 0.) {
						d[i + 1] -= p;
						e[m] = 0.;
						break;
					}
					s = f / r;
					c = g / r;
					g = d[i + 1] - p;
					r = (d[i] - g) * s + 2. * c * b;
					d[i + 1] = g + (p = s * r);
					g = c * r - b;
					if (ifv) {
						f = z[i + 1];
						z[i + 1] = s * z[i] + c * f;
						z[i] = c * z[i] - s * f;
					}
				}
				if (r == 0. && i >= l) continue;
				d[l] -= p;
				e[l] = g;
				e[m] = 0.;
			}
		} while (m != l);
	}
	if (ifv)
		slctsort(d, z);
	else
		sort(d);
	return ierr;
}
