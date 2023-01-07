/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2017-08
*/

#ifndef _PERCENT_H_
#define _PERCENT_H_

#include "stdpfx.h"
#include "iofmt.h"

class Percent {
	friend inline std::ostream &operator<<(std::ostream &os, const Percent &p);
private:
	Real percent;
public:
	Percent(Real ratio) : percent(ratio * 100) {}
	Real operator()() const { return percent; }
};

inline std::ostream &operator<<(std::ostream &os, const Percent &p)
{
	OFS ofs;    ofs.copyfmt(os);
	os << iofmt("fix", 2) << std::setw(5) << p() << "%";
	os.copyfmt(ofs);
	return os;
}

#endif /* _PERCENT_H_ */
