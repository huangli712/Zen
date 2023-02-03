#pragma once

/*
coded by Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 20201106
*/

#include "specs.h"
#include "prmtr.h"

class Green {
public:
	Int nomgs;				// number of positive imaginary frequencies
	Int norbs;				// number of orbitals
    VecCmplx z_omg;
    Vec<MatCmplx> g; 		// green function
  private:
	void read(const Str& file);
public:
	Green(const Prmtr& p, Int nomgs_i, const VecCmplx& z_omg, const Str& file = empty_str);
	Green(Int norbs_i, const Prmtr& p, Int nomgs_i, const VecCmplx& z_omg, const Str& file = empty_str);
	Green(const Green& a, Int orb0, Int orb1);
	Green(const Prmtr& p, Int orb0, Int orb1, const Green& a);
	//virtual Real omg(Int n) const { ERR("omg is used illegally "); return 0; }
	virtual Real omg(Int n) const = 0;
	virtual Str type_info() const {return STR("NOt define");}
	//Cmplx z(Int n) const { ERR("omg is used illegally "); return 0; }
	// virtual Cmplx z(Int n) const = 0;
	virtual Cmplx z(Int n) const = 0;
	inline Green& operator+=(const Green& rhs);
	inline Green& operator-=(const Green& rhs);
	inline Green& operator*=(const Cmplx& c);
	inline MatCmplx& operator[](Idx i);
	inline const MatCmplx& operator[](Idx i) const;
	void write(const Str& green_name, Int iter_cnt=999) const;
	void write(const Str& green_name, const Str& rowname, Int iter_cnt=999) const;
	void write_zen(const Str& green_name, Int nspin=2, Int iter_cnt=999) const;
	void write_zen(const Str& green_name, const Str& rowname, Int nspin=2, Int iter_cnt=999) const;
	inline void reset(Idx n);
};

class ImGreen : public Green {
public:
	Real unit_omg;			// unit imaginary frequency, omg_n = (2 n + 1) unit_omg, 0.01 or 0.02 suggested
private:
	Real sum(const VecCmplx& gf) const;
public:
	ImGreen(const Prmtr& p, const Str& file = empty_str) :
		Green(p, p.num_omg, p.Im_z, file), unit_omg(p.unit_omg) {}
	ImGreen(Int norbs_i, const Prmtr& p, const Str& file = empty_str) :
		Green(norbs_i, p, p.num_omg, p.Im_z, file), unit_omg(p.unit_omg) {}
	ImGreen(const ImGreen& a, Int orb0, Int orb1) :
		Green(a, orb0, orb1), unit_omg(a.unit_omg) {}
	ImGreen(const Prmtr& p, Int orb0, Int orb1, const ImGreen& a) :
		Green(p, orb0, orb1, a), unit_omg(a.unit_omg) {}
	inline ImGreen inverse() const;
	// relative err with b, omg_rsd > 0 is a omg residual
	Real error(const ImGreen& b, Real omg_rsd = 2) const;
	inline ImGreen deduct_iomgn() const;
	// matsubara frequency omg_n = (2 * n + 1) * unit_omg
	// Real omg(Int n) const { return(2 * n + 1) * unit_omg; }
	Real omg(Int n) const {return imag(z_omg[n]);}
    // matsubara frequency i*omg_n
	Cmplx z(Int n) const { return  z_omg[n]; }
	//calculate <c_i^+ c_j> from <<c_j c_i^+>>_(iw_n)
	MatReal particle_number() const;
	void write_occupation_info() const;
	Str type_info() const {return STR("ImGreen");}
};

class ReGreen : public Green {
public:
	Real freq_low;
	Real dlt_freq;
	Real eta_freq;
private:
public:
	ReGreen(const Prmtr& p, const Str& file = empty_str) :
		Green(p, p.nfreq, p.Re_z, file), freq_low(p.freq_low), dlt_freq(p.dlt_freq), eta_freq(p.eta_freq) {}
	ReGreen(Int norbs_i, const Prmtr& p, const Str& file = empty_str) :
		Green(norbs_i, p, p.nfreq, p.Re_z, file), freq_low(p.freq_low), dlt_freq(p.dlt_freq), eta_freq(p.eta_freq) {}
	ReGreen(const ReGreen& a, Int orb0, Int orb1) :
		Green(a, orb0, orb1), freq_low(a.freq_low), dlt_freq(a.dlt_freq), eta_freq(a.eta_freq) {}
	ReGreen(const Prmtr& p, Int orb0, Int orb1, const ReGreen& a) :
		Green(p, orb0, orb1, a), freq_low(p.freq_low), dlt_freq(p.dlt_freq), eta_freq(p.eta_freq) {}
	inline ReGreen inverse() const;
	// Real frequency freq_n = freq_low + n * dlt_freq
	Real omg(Int n) const { return real(z_omg[n]); }
	Cmplx z(Int n) const { return z_omg[n]; }
	Str type_info() const {return STR("ReGreen");}
	
};



Green& Green::operator+=(const Green& rhs)
{
	for_Idx(i, 0, nomgs) (*this)[i] += rhs[i];
	return *this;
}

Green& Green::operator-=(const Green& rhs)
{
	for_Idx(i, 0, nomgs) (*this)[i] -= rhs[i];
	return *this;
}

Green& Green::operator*=(const Cmplx& c)
{
	for_Idx(i, 0, nomgs) (*this)[i] *= c;
	return *this;
}

inline ImGreen operator+(ImGreen lhs, const ImGreen& rhs) { lhs += rhs; return lhs; }
inline ImGreen operator-(ImGreen lhs, const ImGreen& rhs) { lhs -= rhs; return lhs; }
inline ImGreen operator*(Cmplx t, ImGreen rhs) { rhs *= t; return rhs; }

inline ReGreen operator+(ReGreen lhs, const ReGreen& rhs) { lhs += rhs; return lhs; }
inline ReGreen operator-(ReGreen lhs, const ReGreen& rhs) { lhs -= rhs; return lhs; }
inline ReGreen operator*(Cmplx t, ReGreen rhs) { rhs *= t; return rhs; }



MatCmplx& Green::operator[](Idx i)
{
#ifdef _CHECK_BOUNDS_
	if (i >= nomgs) ERR("Vec subscript out of bounds, subscript = " + STR(i) + ", vector size = " + STR(nomgs));
#endif
	return g[i];
}

const MatCmplx& Green::operator[](Idx i) const
{
#ifdef _CHECK_BOUNDS_
	if (i >= nomgs) ERR("Vec subscript out of bounds, subscript = " + STR(i) + ", vector size = " + STR(nomgs));
#endif
	return g[i];
}

ReGreen ReGreen::inverse() const
{
	ReGreen ginv(*this);
	for_Int(i, 0, g.size()) {
		ginv.g[i] = matinvlu(g[i]);
	}
	return ginv;
}

ImGreen ImGreen::inverse() const
{
	ImGreen ginv(*this);
	for_Int(i, 0, g.size()) {
		ginv.g[i] = matinvlu(g[i]);
	}
	return ginv;
}

void Green::reset(Idx n) {
	nomgs = n;
	g.reset(nomgs, MatCmplx(norbs, norbs, 0.));
}


ImGreen ImGreen::deduct_iomgn() const
{
	ImGreen ginveff(*this);
	for_Int(i, 0, g.size()) {
		for_Int(m, 0, g[i].nrows()) {
			ginveff[i][m][m] -= z(i);
		}
	}
	return ginveff;
}