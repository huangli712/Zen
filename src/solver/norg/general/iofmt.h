/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2013 - 2017
*/

#ifndef _IOFMT_H_
#define _IOFMT_H_

#include "stdpfx.h"

// io formatting

static const int w_sec = 10;							// used as setw(w_sec) for default precision seconds output format
static const int w_Real = int(1 - log10(eps_Real)) + 8;	// used as setw(w_Real)     for default high precision output format
static const int w_Real_def = 14;						// used as setw(w_Real_def) for default      precision output format
static const int p_Real = int(1 - log10(eps_Real));	    // used as setprecision(p_Real)     for default high precision output format
static const int p_Real_def = 6;						// used as setprecision(p_Real_def) for default      precision output format

class iofmt {
	friend inline std::ostream &operator<<(std::ostream &o, const iofmt &iof);
	friend inline std::istream &operator>>(std::istream &i, const iofmt &iof);
private:
	Str fmt;
	int p;			// precision
public:
	iofmt(): fmt("def"), p(0) {}
	iofmt(const Str &s, int p_i = 6): fmt(s), p(p_i) {
		if (s == "sec") p = 3;
	}
private:
	std::istream &fmt_def(std::istream &i) const { i.copyfmt(IFS()); return i; }
	std::ostream &fmt_def(std::ostream &o) const { o.copyfmt(OFS()); return o; }
	std::ostream &fmt_sci    (std::ostream &o) const { return fmt_def(o) << std::uppercase << std::scientific << std::setprecision(p_Real); }
	std::ostream &fmt_sci_def(std::ostream &o) const { return fmt_def(o) << std::uppercase << std::scientific << std::setprecision(p_Real_def); }
	std::ostream &fmt_fix    (std::ostream &o) const { return fmt_def(o) << std::fixed << std::setprecision(p); }
};

inline std::ostream &operator<<(std::ostream &o, const iofmt &iof)
{
	if(iof.fmt == "sci")
		return iof.fmt_sci(o);
	else if (iof.fmt == "sci_def")
		return iof.fmt_sci_def(o);
	else if (iof.fmt == "fix" || iof.fmt == "sec")
		return iof.fmt_fix(o);
	else if (iof.fmt == "def")
		return iof.fmt_def(o);
	else {
		WRN("invalid fmt: " + STR(iof.fmt));
		return o;
	}
}

inline std::istream &operator>>(std::istream &i, const iofmt &iof)
{
	if (iof.fmt == "def")
		return iof.fmt_def(i);
	else {
		WRN("invalid fmt: " + STR(iof.fmt));
		return i;
	}
}

inline IFS &biread(IFS &ifs, CharP p, Idx n)
{
	if (n) {
		Idx c = ifs.read(p, n).gcount();
		if (c != n) WRN("c != n, c = " + STR(c) + ", n = " + STR(n));
	}
	return ifs;
}

inline OFS &biwrite(OFS &ofs, CharP p, Idx n)
{
	if (n) {
		ofs.write(p, n);
		ofs.flush();
	}
	return ofs;
}

#endif /* _IOFMT_H_ */
