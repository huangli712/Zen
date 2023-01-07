/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2020-2022
*/
#pragma once

#include "stdpfx.h"
#include "vec.h"
// #include "mat.h"
#include "mympi.h"
//#include <mkl.h>

// class for sparse matrix to save the space, use MKL for CSR saving methond.
// Using the MPI to concle the data in different process.
// BUT now only support for the Real number.

template<typename T>
class SparseMat
    /* When using this special CSR is juest suit for multiplication now */
{
private:
    const MyMpi& mm;			// parameters


    //Int dim;                // the dim of the matrix
    // Three Array Variation of CSR Format : 
    VEC<T> va;           // For the value, store as seperate vector.
    VEC<Int> colnu;      // For the colum number.
    VecInt roin;            // Fot the row Index (default:zero-based indexing)

public:


    Int length_m;           // The Number of rows of the matrix
    Int length_k;           // The Number of columns of the matrix

public:

    SparseMat(const Int& n_i, const MyMpi& mm_i);

    // for the matrix (m-by-k) sperate storge in different PID;
    SparseMat(const Int& m_i, const Int& k_i, const MyMpi& mm_i);

    SparseMat(const Int& length_r, const Int& length_c, const Vec<Int>& roin_i, const Vec<Int>& colnu_i, const Vec<Real>& val_i, const MyMpi& mm_i);

    
    inline void addelement(const Real &value_i, const Int &colnu_i, const Int & roin_i);
    inline void add_empty_row(const Int& roin_i) { if (roin[roin_i + 1] == -1)roin[roin_i + 1] = roin[roin_i]; }
    inline void shrink_to_fit() { va.shrink_to_fit(); colnu.shrink_to_fit(); }

    inline void addelement_one_based(const Real& value_i, const Int& colnu_i, const Int& roin_i);

    inline void ho(Vec<T>& hp, const Vec<T>& p);

    //inline SparseMat<T> operator+(SparseMat<T>& sparmat)const;
    inline SparseMat<T> operator+=(SparseMat<T>& rhs);

    inline Vec<T> operator*(const Vec<T>& KetVec)const;


    //inline void square()const;

    //inline VecReal test(const Vec<T>& KetVec);
    inline SparseMat& operator=(const SparseMat& rhs);					// copy assignment operator, see the remarks above


    //inline void test(Vec<T>& hp, const Vec<T>& KetVec);

    VEC<T> get_va() { return va; }
    VEC<Int> get_colnu() { return colnu; }
    VecInt get_roin() { return roin; }

    void free(){
        std::vector<T>().swap(va);
        std::vector<Int>().swap(colnu);
        roin.reset();
    }

    bool if_hermitian()const;

};


// sparse-matrix types
typedef SparseMat<Int> SparseMatInt;
typedef SparseMat<Idx> SparseMatIdx;
typedef SparseMat<Real> SparseMatReal;
typedef SparseMat<Cmplx> SparseMatCmplx;

// related functions
template<typename T> inline SparseMat<T> operator+(SparseMat<T> lhs,  SparseMat<T> rhs) { return lhs += rhs; }


// zero-based indexing added.
template<typename T>
void SparseMat<T>::addelement(const Real& value_i, const Int& colnu_i, const Int& roin_i)
{
    va.push_back(value_i);
    colnu.push_back(colnu_i);
#ifdef _ASSERTION_
    if (roin[roin_i] == -1)ERR("There are some problem in this row(s)" + NAV2(roin_i, roin[roin_i]));
#endif
    if (roin[roin_i + 1] == -1)roin[roin_i + 1] = roin[roin_i];
    roin[roin_i + 1] += 1;
}

// one-based indexing added.
template<typename T>
void SparseMat<T>::addelement_one_based(const Real& value_i, const Int& colnu_i, const Int& roin_i)
{
    va.push_back(value_i);
    colnu.push_back(colnu_i + 1);
    if (roin[roin_i + 1] == -1)roin[roin_i + 1] = roin[roin_i];
    roin[roin_i + 1] += 1;
}

//Here the KetVec only support for the Real vector.
template<typename T>
void SparseMat<T>::ho(Vec<T>& hp, const Vec<T>& KetVec)
{
    VecPartition vp(mm.np(), mm.id(), length_k);
#ifdef _CHECK_DIMENSION_MATCH_
    ASSERT_EQ(length_k, KetVec.size());
    ASSERT_EQ(length_m, vp.len());
#endif
    VecReal val(va);
    VecReal ret(vp.len(), 0.);
    VecInt colnu_i(colnu);
    if (roin.size() > 1) sparse_MUL(vp.len(), length_k, roin, colnu_i, val, KetVec, ret);
    hp = mm.Allgatherv(ret, vp);
}

template<typename T>
Vec<T> SparseMat<T>::operator*(const Vec<T>& KetVec) const
{
    VecPartition vp(mm.np(), mm.id(), length_k);
#ifdef _CHECK_DIMENSION_MATCH_
    ASSERT_EQ(length_k, KetVec.size());
    ASSERT_EQ(length_m, vp.len());
#endif
    VecReal val(va);
    VecReal ret(vp.len(), 0.);
    VecInt colnu_i(colnu);
    if (roin.size() > 1) sparse_MUL(vp.len(), length_k, roin, colnu_i, val, KetVec, ret);
    Vec<T>hp(mm.Allgatherv(ret, vp));
    return hp;
}

template<typename T>
SparseMat<T>& SparseMat<T>::operator=(const SparseMat& rhs)
{
    va = rhs.va;
    colnu = rhs.colnu;
    roin.reset(rhs.roin);
    length_m = rhs.length_m;
    length_k = rhs.length_k;
    return *this;
}

template<typename T>
SparseMat<T> SparseMat<T>::operator+=(SparseMat<T>& rhs)
{
    VecPartition vp(mm.np(), mm.id(), length_k);
    VecReal val_l(va);
    VecInt colnu_i_l(colnu), roin_l(roin);
    //DBG(NAV3(colnu_i_l.truncate(1254, 1354), val_l.truncate(1254, 1354), roin_l));
    VecReal val_r(rhs.get_va());
    VecInt colnu_i_r(rhs.get_colnu()), roin_r(rhs.get_roin());
    //DBG(NAV3(colnu_i_r, val_r, roin_r));
    VecReal val_ret(va.size(), 0.);
    VecInt roin_ret(roin.size(), 0.), colnu_ret(colnu.size(), 0.);
    sparse_ADD(vp.len(),length_k, roin_l, colnu_i_l, val_l, vp.len(), length_k, roin_r, colnu_i_r, val_r, roin_ret, colnu_ret, val_ret);
    SparseMat<T> retsparmat(vp.len(), length_k, roin_ret, colnu_ret, val_ret, mm);
    //DBG(NAV3(colnu_ret.truncate(1254, 1354), val_ret.truncate(1254, 1354), roin_ret));
    return retsparmat;

}
//template<typename T>
//void SparseMat<T>::square() const
//{
//    //VecPartition vp(mm.np(), mm.id(), length_k);
//    //#ifdef _CHECK_DIMENSION_MATCH_
//    //ASSERT_EQ(length_k, KetVec.size());
//    //ASSERT_EQ(length_m, vp.len());
//    //#endif
//    VecReal val(va);
//    VecReal ret(vp.len(), 0.);
//    VecInt colnu_i(colnu);
//    VecInt roin_ret, colnu_ret;
//    VecReal val_ret;
//    sparse_MUL(vp.len(),length_k,roin,colnu_i,val, vp.len(), length_k, roin, colnu_i, val, roin_ret, colnu_ret, val_ret);
//    //Vec<T>hp(mm.Allgatherv(ret, vp));
//    return hp;
//}

//template<typename T>
//void SparseMat<T>::test(Vec<T>& hp, const Vec<T>& KetVec)
//{
//#ifdef _CHECK_DIMENSION_MATCH_
//    ASSERT_EQ(length_k, KetVec.size());
//    ASSERT_EQ(length_m, length_k);
//    ASSERT_EQ(hp.size(), KetVec.size());
//#endif
//    VecReal val(va);
//    VecInt colnu_i(colnu);
//    sparse_MUL(length_m, roin, colnu_i, val, KetVec, hp);
//}
//
//template<typename T>
//VecReal SparseMat<T>::test(const Vec<T>& KetVec)
//{
//#ifdef _CHECK_DIMENSION_MATCH_
//    ASSERT_EQ(length_k, KetVec.size());
//    ASSERT_EQ(length_m, length_k);
//#endif
//    VecReal val(va);
//    VecReal ret(length_m, 0.);
//    VecInt colnu_i(colnu);
//    sparse_MUL(length_m, roin, colnu_i, val, KetVec, ret);
//    return ret;
//}

template<typename T>
SparseMat<T>::SparseMat(const Int &n_i, const MyMpi& mm_i) : mm(mm_i), length_m(n_i), roin(n_i + 1, 0.)
{
    roin[0] = 0;
}

template<typename T>
inline SparseMat<T>::SparseMat(const Int& m_i, const Int& k_i, const MyMpi& mm_i) : 
    mm(mm_i), length_m(m_i), length_k(k_i), roin(m_i + 1, -1)
{
    roin[0] = 0;

}

template<typename T>
inline SparseMat<T>::SparseMat(const Int& length_r, const Int& length_c, const Vec<Int>& roin_i, const Vec<Int>& colnu_i, const Vec<Real>& val_i, const MyMpi& mm_i) :
    mm(mm_i), length_m(length_r), length_k(length_c), roin(roin_i)
{
    for_Int(i, 0, val_i.size()) {
        colnu.push_back(std::move(colnu_i[i]));
        va.push_back(std::move(val_i[i]));
    }
}

template<typename T>
bool SparseMat<T>::if_hermitian()const
{
#ifdef _CHECK_DIMENSION_MATCH_
    ASSERT_EQ(length_k, length_m); // on support for one phalanx matrix.
#endif
    // zero-based indexing
    bool check(false);
    for_Int(row_i, 0, roin.size()-1) {
        for (Idx clm_i = roin[row_i]; clm_i < roin[row_i+1]; clm_i++){
            Idx x_i(row_i), y_i(colnu[clm_i]);
            Real v_i(va[clm_i]);
            SWAP(x_i, y_i);
            for (Idx n = roin[x_i]; n < roin[x_i+1]; n++) {
                //DBG(NAV4(colnu[n], va[n], y_i, v_i));
                if (colnu[n] == y_i) if( ABS(va[n] - v_i) < 1.E-10) {
                    check = true;
                    break;
                }
                else check = false;
            }
            if (check == false) {
                SWAP(x_i, y_i);
                ERR("Here is not Hermitian!!" + NAV5(v_i, x_i, y_i,va.size(),colnu.size()));
            }
        }
    }
    return check;
}



