/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)date 2013 - 2017
	Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2020 - 2022
*/

#ifndef _STDSFX_H_
#define _STDSFX_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <ios>

#include <string>
#include <vector>
#include <set>
#include <map>

#include <complex>
#include <limits>
#include <algorithm>
#include <functional>
#include <numeric>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cctype>
#include <ctime>
#ifdef _MSC_VER
#include <io.h>
#include <direct.h>
#include <windows.h>
#undef min
#undef max
#else
#include <unistd.h>
#ifdef __APPLE__
        #include <sys/uio.h>
#else
        #include <sys/io.h>
#endif
#include <sys/wait.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif
#include "str.h"
#include "envvar.h"

// pause
inline void suspend()
{
#ifdef _MSC_VER
	system("pause"); 
#endif
}

// exception handling
#define VENUE (empty_str + "[" + __FUNCTION__ + "][" + __FILE__ + "](" + STR(__LINE__) + ")")
#define ERR(msg) { std::cout << "ERR: " << append_space(msg) << VENUE << std::endl; if (debug) suspend(); exit(1); }
#define WRN(msg) { std::cout << "WRN: " << append_space(msg) << VENUE << std::endl; }
#define PIO(msg) { std::cout << "PIO: " << append_space(msg) << std::endl; }
#define DBG(msg) { std::cout << "DBG: " << append_space(msg) << VENUE << std::endl; if (debug) suspend(); }
#define ASSERT_EQ(a,b) { if ((a) != (b)) ERR(STR(#a, " != ", #b) + ", " + NAV2(a, b)); }
#define ASSERT_NE(a,b) { if ((a) == (b)) ERR(STR(#a, " == ", #b) + ", " + NAV2(a, b)); }
#define ASSERT_LT(a,b) { if ((a) >= (b)) ERR(STR(#a, " >= ", #b) + ", " + NAV2(a, b)); }
#define ASSERT_GT(a,b) { if ((a) <= (b)) ERR(STR(#a, " <= ", #b) + ", " + NAV2(a, b)); }
#define ASSERT_LE(a,b) { if ((a) > (b)) ERR(STR(#a, " > ", #b) + ", " + NAV2(a, b)); }
#define ASSERT_GE(a,b) { if ((a) < (b)) ERR(STR(#a, " < ", #b) + ", " + NAV2(a, b)); }

#define ASSERT_EQ1_PS(a0,b0,ps) \
{ \
	if ((a0) != (b0)) \
		ERR(empty_str + #a0 + " != " + #b0 + "\n" + \
		NAVC2(a0, b0) + ps); \
}
#define ASSERT_EQ2_PS(a0,b0,a1,b1,ps) \
{ \
	if ((a0) != (b0) || (a1) != (b1)) \
		ERR(empty_str + #a0 + " != " + #b0 + " || " + #a1 + " != " + #b1 + "\n" + \
		NAVC2(a0, b0) + NAVC2(a1, b1) + ps); \
}
#define ASSERT_EQ3_PS(a0,b0,a1,b1,a2,b2,ps) \
{ \
	if ((a0) != (b0) || (a1) != (b1) || (a2) != (b2)) \
		ERR(empty_str + #a0 + " != " + #b0 + " || " + #a1 + " != " + #b1 + " || " + #a2 + " != " + #b2 + "\n" + \
		NAVC2(a0, b0) + NAVC2(a1, b1) + NAVC2(a2, b2) + ps); \
}
#define ASSERT_EQ4_PS(a0,b0,a1,b1,a2,b2,a3,b3,ps) \
{ \
	if ((a0) != (b0) || (a1) != (b1) || (a2) != (b2) || (a3) != (b3)) \
		ERR(empty_str + #a0 + " != " + #b0 + " || " + #a1 + " != " + #b1 + " || " + #a2 + " != " + #b2 + " || " + #a3 + " != " + #b3 + "\n" + \
		NAVC2(a0, b0) + NAVC2(a1, b1) + NAVC2(a2, b2) + NAVC2(a3, b3) + ps); \
}
#define ASSERT_EQ1(a0,b0) ASSERT_EQ1_PS(a0,b0,empty_str)
#define ASSERT_EQ2(a0,b0,a1,b1) ASSERT_EQ2_PS(a0,b0,a1,b1,empty_str)
#define ASSERT_EQ3(a0,b0,a1,b1,a2,b2) ASSERT_EQ3_PS(a0,b0,a1,b1,a2,b2,empty_str)
#define ASSERT_EQ4(a0,b0,a1,b1,a2,b2,a3,b3) ASSERT_EQ4_PS(a0,b0,a1,b1,a2,b2,a3,b3,empty_str)

#ifdef _MSC_VER
#else
#define nullptr 0
#endif

// basic type names
#ifdef _MSC_VER
typedef __int64 LLInt;
typedef unsigned __int64 ULLInt;
#else
typedef long long int LLInt;
typedef unsigned long long int ULLInt;
#endif
typedef int Int;
typedef unsigned int UInt;
typedef char Char;
typedef char *CharP;
typedef unsigned char UChar;
typedef double Real;
typedef long double LReal;
typedef std::complex<Real> Cmplx;
typedef bool Bool;
typedef std::size_t Idx;		// type for subscript and size of vector and matrix
typedef std::ptrdiff_t PDiff;	// type for difference of subscripts, sizes, or pointers
typedef std::string Str;
typedef std::map<Int, Int> MapInt;
typedef std::map<Idx, Idx> MapIdx;
typedef std::set<Int> SetInt;
typedef std::set<Idx> SetIdx;
typedef std::pair<Int, Int> PairInt;
typedef std::pair<Idx, Idx> PairIdx;
typedef std::pair<Real, Real> PairReal;
typedef std::ifstream IFS;
typedef std::ofstream OFS;
#define VEC std::vector

// macros, there can not be any space between the macro name and the left parenthesis
#define for_Idx(i, bgn, end) for (Idx i = (bgn); i < (end); ++i)
#define for_Int(i, bgn, end) for (Int i = (bgn); i < (end); ++i)
#define Int_for(i, bgn, end) for (Int i = (end-1); i >= (bgn); --i)

// date and time
std::string present();
// date and time with default format
std::string present_time_default_format();

// numerical constants
static const Int max_Int = std::numeric_limits<Int>::max();					// static const Int max_Int = 2147483647;
static const LLInt max_LLInt = std::numeric_limits<LLInt>::max();			// static const LLInt max_LLInt = 9223372036854775807;
static const Cmplx I(0, 1);		// imaginary unit
static const Real  NaN_Real  = std::numeric_limits<Real>::quiet_NaN();
static const Cmplx NaN_Cmplx = std::numeric_limits<Cmplx>::quiet_NaN();
static const Real  eps_Real = std::numeric_limits<Real>::epsilon();
static const Real  min_Real = std::numeric_limits<Real>::min();
static const Real  max_Real = std::numeric_limits<Real>::max();
static const Real  pi_Real = acos(Real(-1.));
static const Real  golden_Real = (sqrt(Real(5.)) - Real(1.)) / Real(2.);
static const Real  eps_double = std::numeric_limits<double>::epsilon();
static const Real  min_double = std::numeric_limits<double>::min();
static const Real  max_double = std::numeric_limits<double>::max();
static const Real  pi_double = acos(-1.);
static const Real  golden_double = (sqrt(5.) - 1.) / 2.;
static const Real  ln_inv_eps = -std::log(eps_Real);

// system command
#ifdef _MSC_VER
static const Str cmd_copy("copy/y");	// copy a file
static const Str cmd_rmf("del/q/s");	// command for remove a file
static const Str cmd_rmd("rd/q/s");		// command for remove a file
static const Str cmd_mkd("md");			// command for mkdir
static const Str cmd_pwd("echo %cd%");	// command for pwd
#else
static const Str cmd_copy("cp");		// copy a file
static const Str cmd_rmf("rm -rf");		// command for remove a file or a directory
static const Str cmd_rmd("rm -rf");		// command for remove a file or a directory
static const Str cmd_mkd("mkdir");		// command for mkdir
static const Str cmd_pwd("pwd");	// command for pwd
#endif

// directory operations
#ifdef _MSC_VER
inline Int dir_exist(const Str& dir) { return _access(dir.c_str(), 0) != -1; }
inline Int MKDIR(const Str& dir) { return _mkdir(dir.c_str()); }
inline void SLEEP(Int milliseconds) { Sleep(milliseconds); }
#else
inline Int dir_exist(const Str& dir) { return access(dir.c_str(), 0) != -1; }
inline Int MKDIR(const Str& dir) { return mkdir(dir.c_str(), 0777); }
inline void SLEEP(Int milliseconds) { usleep(milliseconds * 1000); }
#endif
// get present directory
#ifdef _MSC_VER
inline Str pwd()
{
	char buffer[4096];
	_getcwd(buffer, 4096);
	return Str(buffer);
}
#else
inline Str pwd()
{
	char buffer[4096];
	getcwd(buffer, 4096);
	return Str(buffer);
}
#endif
// print present directory
inline void print_present_dir()
{
	system(cmd_pwd.c_str());
}

// macro-like inline functions
template<typename T> inline T ABS(const T &a) { return a < 0 ? -a : a; }
inline Real ABS(const Cmplx &a) { return std::abs(a); }
template<typename T> inline T EXP(const T &a) { return std::exp(a); }
inline Real EXP(const Int &a) { return std::exp(Real(a)); };
inline Real EXP(const Idx &a) { return std::exp(Real(a)); };
template<typename T> inline T LOG(const T &a) { return std::log(a); }
inline Real LOG(const Int &a) { return std::log(Real(a)); };
inline Real LOG(const Idx &a) { return std::log(Real(a)); };
inline Real ARG(const Real &a) { return std::arg(a); }
inline Real ARG(const Cmplx &a) { return std::arg(a); }
inline Real ARG_php(const Cmplx &a) { Real x = ARG(a); return x + (x <= -pi_Real / 2 ? 2 * pi_Real : 0.); }
template<typename T> inline T SQR(const T &a) { return a * a; }
template<typename T> inline T SQRT(const T &a) { return std::sqrt(a); }
inline Real SQRT(const Int &a) { return std::sqrt(Real(a)); };
inline Real SQRT(const Idx &a) { return std::sqrt(Real(a)); };
inline Real INV(const Real &a) { return 1. / a; }
inline Cmplx INV(const Cmplx &a) { return 1. / a; }
inline Real INV(const Int &a) { return 1. / Real(a); }
inline Real INV(const Idx &a) { return 1. / Real(a); }
inline Real NRMC(const Real &a) { return 1. / SQRT(a); }
inline Real NRMC(const Int &a) { return 1. / SQRT(a); }
inline Real NRMC(const Idx &a) { return 1. / SQRT(a); }
// round to Int
inline Int Int_ROUND(const Real& a) { return Int(a + 0.5); }
// round to LLInt
inline LLInt LLInt_ROUND(const Real& a) { return LLInt(a + 0.5); }
// ceil a so that it is divisible by d
inline Int CEIL(Int a, Int d) { return (a + d - 1) / d * d; }
// ceil a so that it is divisible by d
inline Real CEIL(Real a, Real d) { return std::ceil(a / d) * d; }
// (-1)^n
inline Int minus_one_power(const Int& n) { return 1 - ((n & 1) << 1); }
// (-1)^(m + n)
inline Int minus_one_power(const Int& m, const Int& n) { return minus_one_power(m ^ n); }
// eps accumulated = SQRT(n) * eps, assuming numbers with a relative error of eps
inline Real eps_acc(Idx n, const Real &eps) { return SQRT(n) * eps; }
// eps accumulated, assuming numbers with a relative error 10000 * eps_Real
inline Real eps_acc_def(Idx n) { return 10000. * eps_acc(n, eps_Real); }
// relative error
template<typename T> inline Real RERR(const T &a, const T &b) { return ABS(a - b) * 2 / (ABS(a) + ABS(b) + 1.E-128); }
template<typename T> inline const T &MAX(const T &a, const T &b) { return b > a ? (b) : (a); }
template<typename T> inline const T &MIN(const T &a, const T &b) { return b < a ? (b) : (a); }
template<typename T> inline void SWAP(T &a, T &b) { T dum = a; a = b; b = dum; }
template<typename T> inline int SIGN(const T &a) { return a >= 0 ? 1 : -1; }
template<typename T> inline int SIGN0(const T &a) { return a > 0 ? 1 : (a < 0 ? -1 : 0); }
template<typename T> inline T SIGN(const T &a, const T &b) { return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a); }
template<typename T> inline T SIGN0(const T &a, const T &b) { return b == 0 ? 0 : (b > 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a)); }
template<typename T> inline T MASK(const T &a, const Int &length) { return (T(1) << length) - T(1); }
template<typename T> inline T mod_by_mask(const T &a, const Int &length) { return a & ((T(1) << length) - T(1)); }
// Kronecker delta function
template<typename T> inline Int DELTA(const T &a, const T &b) { return a == b; }

inline Real cnjg(const Real &r) { return r; }
inline Cmplx cnjg(const Cmplx &c) { return Cmplx(c.real(), -c.imag()); }
inline Real real(const Real &v) { return v; }
inline Real imag(const Real &v) { return 0; }
inline Real real(const Cmplx &v) { return v.real(); }
inline Real imag(const Cmplx &v) { return v.imag(); }
inline Cmplx cmplx(const Real &a, const Real &b) { return Cmplx(a, b); }
inline Cmplx cmplx(const Real &a) { return Cmplx(a, 0.); }
inline Real nsqr(const Real &v) { return v * v; }
inline Real nsqr(const Cmplx &v) { return v.real() * v.real() + v.imag() * v.imag(); }
// find the character string length of an object when it is printed out
template<typename T> inline Int TOSTRLEN(T a) { std::ostringstream oss; oss << a; return (Int)oss.str().length(); }
inline Real seconds_since(clock_t bgn) { return Real(clock() - bgn) / CLOCKS_PER_SEC; }
inline void TIME_BGN(const Str &label, clock_t &t_bgn)
{
	std::cout << STR("TIME < " + label + " bgn " + present()) << std::endl;
	t_bgn = clock();
}
inline void TIME_END(const Str &label, const clock_t &t_bgn)
{
	std::cout << STR("TIME > " + label + " end " + present(), "  elapsed ", seconds_since(t_bgn)) << std::endl;
}
inline void print_progress(std::ostream &os, const Str &lbl, Int i, Int n)
{
	Int period = (n - 1) / 100 + 1;
	if ((i + 1) % period) return;

	const int &w_i = TOSTRLEN(n);
	os << lbl + "  " << std::setw(w_i) << i + 1 << " / " << n << "  " << present() << std::endl;
}

// Compute sqrt(a * a + b * b) without destructive underflow or overflow
inline Real pythag(const Real &a, const Real &b)
{
	Real aa = ABS(a), bb = ABS(b);
	return (aa > bb ? aa * sqrt(1. + SQR(bb / aa)) :
		(bb == 0. ? 0. : bb * sqrt(1. + SQR(aa / bb))));
}

// find sqrt(a * a + b * b) without overflow or destructive underflow
// discarded because of its longer code than pythag
// its relative efficency is unknown
inline Real pythag_discarded(const Real &a, const Real &b)
{
	Real p, r, s, t;
	p = MAX(ABS(a), ABS(b));
	if (p == 0.) return p;

	r = SQR(MIN(ABS(a), ABS(b)) / p);
	while (4. + r != 4.) {
		s = r / (4. + r);
		t = 1. + 2 * s;
		p = t * p;
		r = SQR(s / t) * r;
	}
	return p;
}

#endif /* _STDSFX_H_ */
