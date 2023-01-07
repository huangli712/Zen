/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2016-05-13

    copy from mins_ndim.h in Numerical Recipes in C++
    Multidimensional minimization by the Fletcher-Reeves-Polak-Ribiere method
*/

#ifndef _FRPRMN_H_
#define _FRPRMN_H_

#include "toolbox.h"
#include "min1dim.h"

template <class T>
struct Df1dim {
	const VecReal &p;
	const VecReal &xi;
	Int n;
	T &funcd;
	VecReal xt;
	VecReal dft;
	Df1dim(const VecReal &pp, const VecReal &xii, T &funcdd): p(pp),
		xi(xii), n(pp.size()), funcd(funcdd), xt(n), dft(n) {}
	Real operator()(const Real x)
	{
		for (Int j=0;j<n;j++)
			xt[j]=p[j]+x*xi[j];
		return funcd(xt);
	}
	Real df(const Real x)
	{
		Real df1=0.0;
		funcd.df(xt,dft);
		for (Int j=0;j<n;j++)
			df1 += dft[j]*xi[j];
		return df1;
	}
};
template <class T>
struct Dlinemethod {
	VecReal p;
	VecReal xi;
	T &func;
	Int n;
	const Real xtol;
	Dlinemethod(T &funcc, const Real xtoll=3.0e-8): func(funcc), xtol(xtoll) {}
	Real linmin()
	{
		Real ax,xx,xmin;
		n=p.size();
		Df1dim<T> df1dim(p,xi,func);
		ax=0.0;
		xx=1.0;
		Dbrent dbrent(xtol);
		dbrent.bracket(ax,xx,df1dim);
		xmin=dbrent.minimize(df1dim);
		xi *= xmin;
		p += xi;
		return dbrent.fmin;
	}
};
template <class T>
struct Frprmn : Dlinemethod<T> {
	Int iter;
	Real fret;
	using Dlinemethod<T>::func;
	using Dlinemethod<T>::linmin;
	using Dlinemethod<T>::p;
	using Dlinemethod<T>::xi;
	const Int its_max;
	const Real ftol;
	Int cnvtype;
	Frprmn(T &funcd, Int its_max_i = 992, const Real ftoll=3.0e-8): 
		Dlinemethod<T>(funcd),its_max(its_max_i),ftol(ftoll),cnvtype(0) {}
	VecReal minimize(const VecReal &pp)
	{
		const Real EPS=1.0e-18;
		const Real GTOL=1.0e-8;
		Real gg,dgg;
		Int n=pp.size();
		p=pp;
		VecReal g(n),h(n);
		xi.reset(n);
		Real fp=func(p);
		func.df(p,xi);
		h=g=-xi;
		xi=func.rescale_mrn(h);
		for (Int its=0;its<its_max;its++) {
			iter=its+1;
			fret=linmin();
			if (2.0*ABS(fret-fp) <= ftol*(ABS(fret)+ABS(fp)+EPS)) { cnvtype = 1; return p; }
			fp=fret;
			func.df(p,xi);
			Real test=0.0;
			Real den=MAX(fp,1.0);
			for (Int j=0;j<n;j++) {
				Real temp=ABS(xi[j])*MAX(ABS(p[j]),1.0)/den;
				if (temp > test) test=temp;
			}
			if (test < GTOL) { cnvtype = 2; return p; }
			gg = DOT(g, g);
			dgg=0.0;
			for (Int j=0;j<n;j++) {
//				dgg += xi[j]*xi[j];
				dgg += (xi[j]+g[j])*xi[j];
			}
			if (gg == 0.0)  { cnvtype = 3; return p; }
			g=-xi;
			h=g+(dgg/gg)*h;
			xi=func.rescale_mrn(h);
		}
		return p;
	}
};

#endif /* _FRPRMN_H_ */
