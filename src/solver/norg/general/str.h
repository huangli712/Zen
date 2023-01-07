/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)date 2013 - 2017
	Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2020 - 2022
*/

#ifndef _TOSTR_H_
#define _TOSTR_H_

#include <iostream>
#include <iomanip>
#include <string>

const std::string empty_str;

inline std::string STR() { return std::string(); }
inline std::string append_space(const std::string &s)
{
	if (s.empty()) return s;
	return s[s.size() - 1] == '\n' ? s : s + " ";
}

template<typename T0>                                                                                                                                                                                                    inline std::string STR(const T0 &t0)                                                                                                                                                                                                                   { std::ostringstream ret; ret << t0;                                                                                           return ret.str(); }
template<typename T0, typename T1>                                                                                                                                                                                       inline std::string STR(const T0 &t0, const T1 &t1)                                                                                                                                                                                                     { std::ostringstream ret; ret << t0 << t1;                                                                                     return ret.str(); }
template<typename T0, typename T1, typename T2>                                                                                                                                                                          inline std::string STR(const T0 &t0, const T1 &t1, const T2 &t2)                                                                                                                                                                                       { std::ostringstream ret; ret << t0 << t1 << t2;                                                                               return ret.str(); }
template<typename T0, typename T1, typename T2, typename T3>                                                                                                                                                             inline std::string STR(const T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3)                                                                                                                                                                         { std::ostringstream ret; ret << t0 << t1 << t2 << t3;                                                                         return ret.str(); }
template<typename T0, typename T1, typename T2, typename T3, typename T4>                                                                                                                                                inline std::string STR(const T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4)                                                                                                                                                           { std::ostringstream ret; ret << t0 << t1 << t2 << t3 << t4;                                                                   return ret.str(); }
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>                                                                                                                                   inline std::string STR(const T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5)                                                                                                                                             { std::ostringstream ret; ret << t0 << t1 << t2 << t3 << t4 << t5;                                                             return ret.str(); }
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>                                                                                                                      inline std::string STR(const T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6)                                                                                                                               { std::ostringstream ret; ret << t0 << t1 << t2 << t3 << t4 << t5 << t6;                                                       return ret.str(); }
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>                                                                                                         inline std::string STR(const T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6, const T7 &t7)                                                                                                                 { std::ostringstream ret; ret << t0 << t1 << t2 << t3 << t4 << t5 << t6 << t7;                                                 return ret.str(); }
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>                                                                                            inline std::string STR(const T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8)                                                                                                   { std::ostringstream ret; ret << t0 << t1 << t2 << t3 << t4 << t5 << t6 << t7 << t8;                                           return ret.str(); }
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>                                                                               inline std::string STR(const T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8, const T9 &t9)                                                                                     { std::ostringstream ret; ret << t0 << t1 << t2 << t3 << t4 << t5 << t6 << t7 << t8 << t9;                                     return ret.str(); }
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename Ta>                                                                  inline std::string STR(const T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8, const T9 &t9, const Ta &ta)                                                                       { std::ostringstream ret; ret << t0 << t1 << t2 << t3 << t4 << t5 << t6 << t7 << t8 << t9 << ta;                               return ret.str(); }
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename Ta, typename Tb>                                                     inline std::string STR(const T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8, const T9 &t9, const Ta &ta, const Tb &tb)                                                         { std::ostringstream ret; ret << t0 << t1 << t2 << t3 << t4 << t5 << t6 << t7 << t8 << t9 << ta << tb;                         return ret.str(); }
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename Ta, typename Tb, typename Tc>                                        inline std::string STR(const T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8, const T9 &t9, const Ta &ta, const Tb &tb, const Tc &tc)                                           { std::ostringstream ret; ret << t0 << t1 << t2 << t3 << t4 << t5 << t6 << t7 << t8 << t9 << ta << tb << tc;                   return ret.str(); }
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename Ta, typename Tb, typename Tc, typename Td>                           inline std::string STR(const T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8, const T9 &t9, const Ta &ta, const Tb &tb, const Tc &tc, const Td &td)                             { std::ostringstream ret; ret << t0 << t1 << t2 << t3 << t4 << t5 << t6 << t7 << t8 << t9 << ta << tb << tc << td;             return ret.str(); }
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename Ta, typename Tb, typename Tc, typename Td, typename Te>              inline std::string STR(const T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8, const T9 &t9, const Ta &ta, const Tb &tb, const Tc &tc, const Td &td, const Te &te)               { std::ostringstream ret; ret << t0 << t1 << t2 << t3 << t4 << t5 << t6 << t7 << t8 << t9 << ta << tb << tc << td << te;       return ret.str(); }
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename Ta, typename Tb, typename Tc, typename Td, typename Te, typename Tf> inline std::string STR(const T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8, const T9 &t9, const Ta &ta, const Tb &tb, const Tc &tc, const Td &td, const Te &te, const Tf &tf) { std::ostringstream ret; ret << t0 << t1 << t2 << t3 << t4 << t5 << t6 << t7 << t8 << t9 << ta << tb << tc << td << te << tf; return ret.str(); }


#define NAME(a) STR(#a)

#define NAV(a) STR(#a, " = ", a)

#define NAV1(a0)                         NAV(a0)
#define NAV2(a0,a1)                      NAV(a0) + ", " + NAV(a1)
#define NAV3(a0,a1,a2)                   NAV(a0) + ", " + NAV(a1) + ", " + NAV(a2)
#define NAV4(a0,a1,a2,a3)                NAV(a0) + ", " + NAV(a1) + ", " + NAV(a2) + ", " + NAV(a3)
#define NAV5(a0,a1,a2,a3,a4)             NAV(a0) + ", " + NAV(a1) + ", " + NAV(a2) + ", " + NAV(a3) + ", " + NAV(a4)
#define NAV6(a0,a1,a2,a3,a4,a5)          NAV(a0) + ", " + NAV(a1) + ", " + NAV(a2) + ", " + NAV(a3) + ", " + NAV(a4) + ", " + NAV(a5)
#define NAV7(a0,a1,a2,a3,a4,a5,a6)       NAV(a0) + ", " + NAV(a1) + ", " + NAV(a2) + ", " + NAV(a3) + ", " + NAV(a4) + ", " + NAV(a5) + ", " + NAV(a6)
#define NAV8(a0,a1,a2,a3,a4,a5,a6,a7)    NAV(a0) + ", " + NAV(a1) + ", " + NAV(a2) + ", " + NAV(a3) + ", " + NAV(a4) + ", " + NAV(a5) + ", " + NAV(a6) + ", " + NAV(a7)
#define NAV9(a0,a1,a2,a3,a4,a5,a6,a7,a8) NAV(a0) + ", " + NAV(a1) + ", " + NAV(a2) + ", " + NAV(a3) + ", " + NAV(a4) + ", " + NAV(a5) + ", " + NAV(a6) + ", " + NAV(a7) + ", " + NAV(a8)

#define NAVC1(a0)                         NAV(a0) + "\n"
#define NAVC2(a0,a1)                      NAV(a0) + "\n" + NAV(a1) + "\n"
#define NAVC3(a0,a1,a2)                   NAV(a0) + "\n" + NAV(a1) + "\n" + NAV(a2) + "\n"
#define NAVC4(a0,a1,a2,a3)                NAV(a0) + "\n" + NAV(a1) + "\n" + NAV(a2) + "\n" + NAV(a3) + "\n"
#define NAVC5(a0,a1,a2,a3,a4)             NAV(a0) + "\n" + NAV(a1) + "\n" + NAV(a2) + "\n" + NAV(a3) + "\n" + NAV(a4) + "\n"
#define NAVC6(a0,a1,a2,a3,a4,a5)          NAV(a0) + "\n" + NAV(a1) + "\n" + NAV(a2) + "\n" + NAV(a3) + "\n" + NAV(a4) + "\n" + NAV(a5) + "\n"
#define NAVC7(a0,a1,a2,a3,a4,a5,a6)       NAV(a0) + "\n" + NAV(a1) + "\n" + NAV(a2) + "\n" + NAV(a3) + "\n" + NAV(a4) + "\n" + NAV(a5) + "\n" + NAV(a6) + "\n"
#define NAVC8(a0,a1,a2,a3,a4,a5,a6,a7)    NAV(a0) + "\n" + NAV(a1) + "\n" + NAV(a2) + "\n" + NAV(a3) + "\n" + NAV(a4) + "\n" + NAV(a5) + "\n" + NAV(a6) + "\n" + NAV(a7) + "\n"
#define NAVC9(a0,a1,a2,a3,a4,a5,a6,a7,a8) NAV(a0) + "\n" + NAV(a1) + "\n" + NAV(a2) + "\n" + NAV(a3) + "\n" + NAV(a4) + "\n" + NAV(a5) + "\n" + NAV(a6) + "\n" + NAV(a7) + "\n" + NAV(a8) + "\n"

template<typename T> inline std::string left_justify(T a, const int &width)
{
	std::string a_str = STR(a);
	int space_len = width - (int)a_str.size();
	space_len = space_len < 0 ? 0 : space_len;
	return a_str + std::string(space_len, ' ');
}
template<typename T> inline std::string rght_justify(T a, const int &width, const char filled_with = ' ')
{
	std::string a_str = STR(a);
	int space_len = width - (int)a_str.size();
	space_len = space_len < 0 ? 0 : space_len;
	return std::string(space_len, filled_with) + a_str;
}
template<typename T> inline std::string prefill0(T a, const int& width)
{
	std::string a_str = STR(a);
	int space_len = width - (int)a_str.size();
	space_len = space_len < 0 ? 0 : space_len;
	return std::string(space_len, '0') + a_str;
}
#endif /* _TOSTR_H_ */
