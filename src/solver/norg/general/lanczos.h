#pragma once

/*
code developed by
    Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2020 - 2022
*/


#include "miblas.h"
#include "milapacke.h"		// to get the trd_heevr_vd(diagonalize Tridiagonal matrix function)
#include "stdpfx.h"
#include "vec.h"
#include "mat.h"
#include "random.h"

template<typename T>
bool residual_method_for_b(VEC<T> ltd, VEC<T> lt_sd, Real b, Int k, bool fast_mode = false)
{
    using namespace std;
    Mat<T> ev(ltd.size(), ltd.size(), 0.);
    Vec<T> va(ltd);
    Vec<T> vb(lt_sd);
    Int info = trd_heevr_qr(va, vb, ev); if (info != 0)ERR("the i-th parameter had an illegal value." + NAV3(info, va, vb));
    if ((fast_mode) && ABS(b * ev[k][0]) < 1.E-10) return true;
    if (!(fast_mode) && ABS(b * ev[k][0]) < 1.E-12) return true;
    return false;
}
template<typename T>
inline Real compare_error(T& pre, T& lst) { return (2 * ABS(lst - pre)) / (ABS(pre) + ABS(lst) + 1.); }

template<typename T>
void orthogonalization_one_to_multiple(Vec<T>& v, const VEC<Vec<T>>& m, Int e)
//	orthogonalization v to m[0] through m[e - 1], which are assumed to be orthonormal
{
    Real mag0, mag1;
    mag1 = sqrt(real(DOT(v, v)));		// mag1 != 0 is supposed
    mag0 = mag1 * 1.E+16;
    while (mag1 / mag0 < 0.98) {
        for_Int(i, 0, e) v -= DOT(m[i], v) * m[i];
        mag0 = mag1;
        mag1 = sqrt(real(DOT(v, v)));		// mag1 != 0 is supposed
    }
}
//template<typename T>
//inline void orthonormalize_one_to_multiple(Vec<T>& v, const VEC<Vec<T>>& m, Int e){ for_Idx(i, 0, e) v -= DOT(m[i], v) * m[i]; }

template<typename T, typename F>
VecInt lanczos(VecReal& evals, Mat<T>& evecs, VecInt& ev_dgcy, Idx n, Int wish_nev, F& H,
    VecReal& inital_state, bool fast_mode = false, Int max_krylov = 9999)
    //	Lanczos algorithm to find several lowest eigenpairs of a Hermitian matrix
    //  only typename F(matrix) H needed with well define operator*(which to returns H * |p>)
    //	typename T can only be Real or Cmplx(NOT USED/TESTED)
    //Output:
    //  eval                is the eigenvalues                              (Dynamic size)
    //  evec                is the eigenvectors, one row contains one vector(Dynamic nrow)
    //  ev_dgcy             a set of energy level's degeneracy              (size fixed)
    //Input:
    //  n                   is the dimension of the H's space.
    //  wish_nev            is the expected maximum eigenstates to be found (suggest: 1).
    //  H                   typename F(matrix) H needed with well define operator*(which to returns H * |p>)
    //  inital_state        (not necessary) if set it as blank, Then we will generate a random vector.
    //  max_krylov          (not necessary) is the maximum size of the Krylov space

{
    // if (fast_mode) PIO("Dear customer, you are now in the fast mode, in this mode you only can find the ONE eigen pair.");
#ifdef _CHECK_DIMENSION_MATCH_
    //ASSERT_EQ(n, evec[].());
    //ASSERT_EQ(eval.size(), wish_nev);
    //ASSERT_EQ(evec.size(), wish_nev);
    ASSERT_EQ(ev_dgcy.size(), wish_nev);
#endif
    using namespace std;
    Int nev(wish_nev), e(0);
    VEC<VecReal> evec;
    VEC<Real> eval;
    VEC<Int> Krylovsize;
    while (true)
    {
        // for tridiagonal matrix
        VEC<Real> ltd;	    // diagonal elements / eigenvalues
        VEC<Real> lt_sd;	// sub-diagonal elements
        // prepare an initial Lanczos vector
        VecReal temp = inital_state.normalize();
        if (e > 0) orthogonalization_one_to_multiple(temp, evec, e);
        // for_Idx(i, 0, e) v0 -= DOT(evec[i], v0) * evec[i];
        temp.normalize();
        VecReal v0 = temp;
        VecReal v1(H * v0);
        Idx k = 0;
        Real a(0.), b(0.);
        while (true)
        {
            a = DOT(v1, v0);
            ltd.push_back(a);
            if (!(fast_mode)) if (k >= n/10 || k >= 60 ) if (residual_method_for_b(ltd, lt_sd, b, k, fast_mode)) break;
            if ((fast_mode))  if (k >= n/10 || k >= 60 ) if (residual_method_for_b(ltd, lt_sd, b, k, fast_mode)) break;
            if (k >= max_krylov) {
                WRN("Getting Wrong with someting? krylov was too large!" + NAV2(e, k));
                if (k >= 20 + max_krylov) break;
            }
            k++;
            v1 -= a * v0;
            b = v1.norm();
            lt_sd.push_back(b);
            if (e > 0) orthogonalization_one_to_multiple(v1, evec, e);
            //for_Idx(i, 0, e) v1 -= DOT(evec[i], v1) * evec[i];
            v1.normalize(); // v1 *= (1 / b);
            SWAP(v0, v1);
            v1 *= -b;
            v1 += H * v0;
        }
        // WRN(NAV(k) + "The iteration" + NAV(e));
        VecReal va(ltd), vb(lt_sd);
#ifdef _ASSERTION_
        if (va.size() - 1 != vb.size())ERR("Wrong with lanczos" + NAV3(e, va.size(), vb.size()));
#endif
        MatReal ev(ltd.size(), ltd.size(), 0.);
        trd_heevr_qr(va, vb, ev);
        VecReal gs_kvsp(ev.tr()[0]);                //ground state at Krylov space.
        VecReal gs(inital_state.size(), 0.);
        v0 = temp;
        gs += gs_kvsp[0] * v0;
        v1 = H * v0;
        k = 1;
        while (true)
        {
            v1 -= ltd[k - 1] * v0;
            if (e > 0) orthogonalization_one_to_multiple(v1, evec, e);
            // for_Idx(i, 0, e) v1 -= DOT(evec[i], v1) * evec[i];
            if (k == gs_kvsp.size() - 1) break;
            v1.normalize();
            gs += gs_kvsp[k] * v1;
            SWAP(v0, v1);
            v1 *= -lt_sd[k - 1];
            v1 += H * v0;  k++;
        }
        Krylovsize.push_back(k);
        if (e > 0 && compare_error(eval[e - 1], va[0]) < 1.E-8) nev++;
            // eval.reset(nev, eval.begin());
            // eval.push_back()
            // if(fast_mode) WRN("The eigen value may have degeneracy. This may take some attention.");
        //WRN(NAV5(wish_nev,e, nev,eval.size(),evec.size()));
        // std::cout << "The eigenvalue" << iofmt("sci") << va[0] << std::endl;
        if (e >= nev) break;
        Idx level = wish_nev + e - nev;
        evec.push_back(gs);
        eval.push_back(va[0]); ev_dgcy[level]++;
        e++;
        if (fast_mode) break;
    }
    evecs.reset(evec);
    evals.reset(eval);
    VecInt krylov_size(Krylovsize);
    return krylov_size;
}

//template<typename T, typename F>
//void lanczos(VecReal& eval, Mat<T>& evec, VecInt& ev_dgcy, Int wish_nev, F& H, Int max_krylov = 999)
//{
//    UURandomSFMT uur;
//    Vec <T> ev(evec[e]);
//    if (giv <= e) { uur(ev); ev -= 0.5;}
//}

