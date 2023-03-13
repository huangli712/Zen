#pragma once
/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2023.03.03
*/


#include <bitset>
#include <string>
#include"norg.h"

// adaptive sampling natural configuration interaction (Asnci)
// natural configuration interaction (Nci)

// typedef pair<VEC<Str>,VEC<Real>> Nci;
// typedef pair<VEC<__uint128_t>,VEC<Real>> Nci;
typedef pair<VEC<UInt*>,VEC<Real>> Nci;
class Asnci 
{
	const MyMpi& mm;				// parameters
	const Prmtr& p;					// parameters
	const NocSpace& nosp;			// parameters

    // using namespace std;
    // VEC<Str>            cfig_inital;
    // VEC<Real>           rank_inital;
    // Vec<Str>            cfig_trncat;
    // Vec<Real>           rank_trncat;
    // Nci inital;    

    Real groundE;
    Nci trncat;
    std::map<UInt*, Idx>  cfig_idx;

    MatReal hop_h;
    VEC<Int> mayhop;

public:
    Idx dim;                // The truncated space size

private:

    VEC<Int> find_mayhop();
    
    Nci git_nci(const NORG& norg);
    
    void expand(Nci& natural_cfgs);

    Str to_binary_string(unsigned long num) {
        using namespace std;
        bitset<sizeof(unsigned long) * 8> binary(num);
        return binary.to_string();
    }

    bool judge(Str& cfg_str) {
        for_Idx(i, 0, cfg_str.size()) if((cfg_str[i] != '1' && cfg_str[i] != '0')) return false;
        return true;
    }

    bool judge(const Str& cfg_str,const Int& pos) {
        bool flag(true);
        Int crt(pos/hop_h.ncols()), ann(pos%hop_h.ncols());
        flag = (cfg_str[ann]=='1') ? true : false;
        flag = (cfg_str[crt]=='0') ? true : false;
        return flag;
    }

    Str change_cfg_str(const Str& cfg_str, Int& pos) {
        Int crt(pos/hop_h.ncols()), ann(pos%hop_h.ncols());
        Str temp_c = cfg_str;
        temp_c[ann] = (temp_c[ann]=='1') ? '0' : 'x';
        temp_c[crt] = (temp_c[crt]=='0') ? '1' : 'x';
        return temp_c;
    }

    // The hamilton form:
    Real hamilton_value(Str alpha, Str beta_i = Str()) {
        Real value(0.);
        if(beta_i.length() == 0) {
            value += 
            return value;
        }
        else {
            value += 
            return value;
        }
    }

    Real cfi2rank(Str alpha, Vec<Str> beta) {
        Real rank(0.);
        for_Idx(i, 0, beta.size()){
            Real upper = hamilton_value(alpha, beta[i]);
            rank += upper / (groundE - hamilton_value(alpha));
        }
        return rank;
    }

    Nci truncation(const Nci& inital);

    Tab find_table(Str inter_type);

public:

// Asnci: mode = 0: assume NO converged;  
//        mode = 1: NO converg with adaptive sampling configuration interaction.
Asnci(const NORG& norg, Idx trncat_size, const Int mode = 0);

// out put the NORG class with the table
NORG get_norg(Tab table, Int mode = 0);


};