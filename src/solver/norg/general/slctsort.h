/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2013 - 2017

	selection sort
	it is easy, but slow. don't use it for large-scale or even moderate-scale problems.
*/
#ifndef _SLCTSORT_H_
#define _SLCTSORT_H_

#include "stdpfx.h"
#include "vec.h"
#include "mat.h"

#define SLCTSORT(Td,d,z,zn,cmp) \
{ \
	ASSERT_EQ(d.size(), zn); \
	Idx n = d.size(); \
	for_Idx (i, 0, n) { \
		Td p = d[i]; \
		Idx k = i; \
		for_Idx (j, i + 1, n) if (cmp(d[j], p)) p = d[k = j]; \
		if (k != i) { \
			SWAP(d[i], d[k]); \
			SWAP(z[i], z[k]); \
		} \
	} \
}

#define SLCTSORT2(Td,d,z0,z0n,z1,z1n,cmp) \
{ \
	ASSERT_EQ(d.size(), z0n); \
	ASSERT_EQ(d.size(), z1n); \
	Idx n = d.size(); \
	for_Idx (i, 0, n) { \
		Td p = d[i]; \
		Idx k = i; \
		for_Idx (j, i + 1, n) if (cmp(d[j], p)) p = d[k = j]; \
		if (k != i) { \
			SWAP(d[i], d[k]); \
			SWAP(z0[i], z0[k]); \
			SWAP(z1[i], z1[k]); \
		} \
	} \
}

/* ------------------------------------------ base macro end ------------------------------------------ */

#ifdef _CPP11_

template<typename Td, typename F = std::less<Td> > inline void sort(VEC<Td> &d, const F &cmp = F()) { sort(d.begin(), d.end(), cmp); }
template<typename Td, typename F = std::less<Td> > inline void sort(Vec<Td> &d, const F &cmp = F()) { sort(d.begin(), d.end(), cmp); }

#define TEMPLATE_SLCTSORT(Vd,Vz,z,zn) \
template<typename Td, typename Tz, typename F = std::less<Td> > \
inline void slctsort(Vd<Td> &d, Vz<Tz> &z, const F &cmp = F()) \
{ SLCTSORT(Td, d, z, zn, cmp); }

TEMPLATE_SLCTSORT(VEC, VEC, z, z.size());
TEMPLATE_SLCTSORT(VEC, Vec, z, z.size());
TEMPLATE_SLCTSORT(Vec, VEC, z, z.size());
TEMPLATE_SLCTSORT(Vec, Vec, z, z.size());
TEMPLATE_SLCTSORT(VEC, Mat, z, z.nrows());
TEMPLATE_SLCTSORT(Vec, Mat, z, z.nrows());

template<typename Td, typename Tz, typename F = std::less<Td> >
inline void slctsort_matcol(Vec<Td> &d, Mat<Tz> &z, const F &cmp = F())
{
	Mat<Tz> zt = z.tr();
	SLCTSORT(Td, d, zt, z.ncols(), cmp);
	z = zt.tr();
}

#define TEMPLATE_SLCTSORT2(Vd,Vz0,z0,z0n,Vz1,z1,z1n) \
template<typename Td, typename Tz0, typename Tz1, typename F = std::less<Td> > \
inline void slctsort(Vd<Td> &d, Vz0<Tz0> &z0, Vz1<Tz1> &z1, const F &cmp = F()) \
{ SLCTSORT2(Td, d, z0, z0n, z1, z1n, cmp); }

TEMPLATE_SLCTSORT2(VEC, VEC, z0, z0.size(), VEC, z1, z1.size());

/* ------------------------------------------ C++11 end ------------------------------------------ */

#else

template<typename Td> inline void sort(VEC<Td> &d) { std::less<Td> cmp; sort(d.begin(), d.end(), cmp); }
template<typename Td> inline void sort(Vec<Td> &d) { std::less<Td> cmp; sort(d.begin(), d.end(), cmp); }

#define TEMPLATE_SLCTSORT(Vd,Vz,z,zn) \
template<typename Td, typename Tz> \
inline void slctsort(Vd<Td> &d, Vz<Tz> &z) \
{ std::less<Td> cmp; SLCTSORT(Td, d, z, zn, cmp); }

// use std::vector d as a reference to sort std::vector d and std::vector z
TEMPLATE_SLCTSORT(VEC, VEC, z, z.size());
// use std::vector d as a reference to sort std::vector d and Vec z
TEMPLATE_SLCTSORT(VEC, Vec, z, z.size());
// use Vec d as a reference to sort Vec d and std::vector z
TEMPLATE_SLCTSORT(Vec, VEC, z, z.size());
// use Vec d as a reference to sort Vec d and Vec z
TEMPLATE_SLCTSORT(Vec, Vec, z, z.size());

// use std::vector d as a reference to sort std::vector d and Mat z
// the rows of z are treated as elements
TEMPLATE_SLCTSORT(VEC, Mat, z, z.nrows());

// use Vec d as a reference to sort Vec d and Mat z
// the rows of z are treated as elements
TEMPLATE_SLCTSORT(Vec, Mat, z, z.nrows());



// use Vec d as a reference to sort Vec d and Mat z
// the columns of z are treated as elements
template<typename Td, typename Tz>
inline void slctsort_matcol(Vec<Td> &d, Mat<Tz> &z)
{
	Mat<Tz> zt = z.tr();
	std::less<Td> cmp;
	SLCTSORT(Td, d, zt, z.ncols(), cmp);
	z = zt.tr();
}

#define TEMPLATE_SLCTSORT2(Vd,Vz0,z0,z0n,Vz1,z1,z1n) \
template<typename Td, typename Tz0, typename Tz1> \
inline void slctsort(Vd<Td> &d, Vz0<Tz0> &z0, Vz1<Tz1> &z1) \
{ std::less<Td> cmp; SLCTSORT2(Td, d, z0, z0n, z1, z1n, cmp); }

// use std::vector d as a reference to sort std::vector d, z0, and z1
TEMPLATE_SLCTSORT2(VEC, VEC, z0, z0.size(), VEC, z1, z1.size());

/* ------------------------------------------ cmp ------------------------------------------ */

template<typename Td, typename F> inline void sort(VEC<Td> &d, const F &cmp) { sort(d.begin(), d.end(), cmp); }
template<typename Td, typename F> inline void sort(Vec<Td> &d, const F &cmp) { sort(d.begin(), d.end(), cmp); }

#define TEMPLATE_SLCTSORT_CMP(Vd,Vz,z,zn) \
template<typename Td, typename Tz, typename F> \
inline void slctsort(Vd<Td> &d, Vz<Tz> &z, const F &cmp) \
{ SLCTSORT(Td, d, z, zn, cmp); }

// use std::vector d as a reference to sort std::vector d and std::vector z with comparer cmp
TEMPLATE_SLCTSORT_CMP(VEC, VEC, z, z.size());
// use std::vector d as a reference to sort std::vector d and Vec z with comparer cmp
TEMPLATE_SLCTSORT_CMP(VEC, Vec, z, z.size());
// use Vec d as a reference to sort Vec d and std::vector z with comparer cmp
TEMPLATE_SLCTSORT_CMP(Vec, VEC, z, z.size());
// use Vec d as a reference to sort Vec d and Vec z with comparer cmp
TEMPLATE_SLCTSORT_CMP(Vec, Vec, z, z.size());

// use std::vector d as a reference to sort std::vector d and Mat z with comparer cmp
// the rows of z are treated as elements
TEMPLATE_SLCTSORT_CMP(VEC, Mat, z, z.nrows());

// use Vec d as a reference to sort Vec d and Mat z with comparer cmp
// the rows of z are treated as elements
TEMPLATE_SLCTSORT_CMP(Vec, Mat, z, z.nrows());



// use Vec d as a reference to sort Vec d and Mat z with comparer cmp
// the columns of z are treated as elements
template<typename Td, typename Tz, typename F>
inline void slctsort_matcol(Vec<Td> &d, Mat<Tz> &z, const F &cmp)
{
	Mat<Tz> zt = z.tr();
	SLCTSORT(Td, d, zt, z.ncols(), cmp);
	z = zt.tr();
}



#define TEMPLATE_SLCTSORT2_CMP(Vd,Vz0,z0,z0n,Vz1,z1,z1n) \
template<typename Td, typename Tz0, typename Tz1, typename F> \
inline void slctsort(Vd<Td> &d, Vz0<Tz0> &z0, Vz1<Tz1> &z1, const F &cmp) \
{ SLCTSORT2(Td, d, z0, z0n, z1, z1n, cmp); }

// use std::vector d as a reference to sort std::vector d, z0, and z1 with comparer cmp
TEMPLATE_SLCTSORT2_CMP(VEC, VEC, z0, z0.size(), VEC, z1, z1.size());

#endif

#endif /* _SLCTSORT_H_ */
