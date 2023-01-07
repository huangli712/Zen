#include "comdivs.h"

/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2021-2022
*/

const BnmlFast ComDivs::bf = BnmlFast();

using namespace std;

ComDivs::ComDivs(Int i_i, VecInt ne_i, VecInt ns_i, bool quick) : idx(i_i), ne(ne_i), ns(ns_i), cf(cf_ComDivs(i_i))
{
	if (!quick)
	{
		div_orb_e = std::get<0>(find_div_orb());
		div_orb_c = std::get<1>(find_div_orb());
	}
}

ComDivs::ComDivs(Int i_i, MatInt ne_i, MatInt ns_i, bool quick) : idx(i_i), ne(ne_i.vec()), ns(ns_i.vec()), cf(cf_ComDivs(i_i))
{
	if (!quick)
	{
		div_orb_e = std::get<0>(find_div_orb());
		div_orb_c = std::get<1>(find_div_orb());
	}
}

ComDivs::ComDivs(Int i_i, VecInt ne_i, VecInt ns_i) : idx(i_i), ne(ne_i), ns(ns_i), cf(cf_ComDivs(i_i)), div_orb_e(std::get<0>(find_div_orb())), div_orb_c(std::get<1>(find_div_orb()))
{}

ComDivs::ComDivs(Int i_i, MatInt ne_i, MatInt ns_i) : idx(i_i), ne(ne_i.vec()), ns(ns_i.vec()), cf(cf_ComDivs(i_i)), div_orb_e(std::get<0>(find_div_orb())), div_orb_c(std::get<1>(find_div_orb()))
{}

ComDivs::ComDivs(VecOnb &newcf, VecInt ne_i, VecInt ns_i) : ne(ne_i), ns(ns_i), cf(newcf)
{
	if (ne_i.size() != ns_i.size())
		ERR("some err with input" + NAV2(ne_i, ns_i));
	idx = find_cf_idx(cf);
}

ComDivs::ComDivs(VecOnb &newcf, MatInt ne_i, MatInt ns_i) : ne(ne_i.vec()), ns(ns_i.vec()), cf(newcf)
{
	if (ne_i.size() != ns_i.size())
		ERR("some err with input" + NAV2(ne_i, ns_i));
	idx = find_cf_idx(cf);
}

VecOnb ComDivs::cf_ComDivs(const Int idx_in)
{
	VecInt base(ns.size() + 1, 1);
	for_Int(i, 1, base.size()) base[i] = base[i - 1] * bf(ns[i - 1], ne[i - 1]);
	VecInt rep(FreeBaseDeCode(idx_in, base));
	VecOnb  a(ns.size());
	for_Int(i, 0, a.size()) {
		Onb a_i(ne[i], ns[i], rep[i]);
		a[i] = a_i;
	}
	return a;
}

ComDivs::OrbTuple ComDivs::find_div_orb()
{
	OrbTuple div_orb;
	vector<VecInt> div_orb_es, div_orb_cs;
	for_Int(i, 0, ns.size()) {
		Int k_e(0), k_c(0);
		VecInt div_orb_e(ne[i], -1), div_orb_c((ns[i] - ne[i]), -1);
		for_Int(j, 0, ns[i]) {
			if (cf[i].isocc(j)) div_orb_e[k_e++] = j;
			else if (cf[i].isuno(j)) div_orb_c[k_c++] = j;
			else ERR("Wrong  with count orbit idx: ");
		}
#ifdef _CHECK_BOUNDS_
		for_Int(i, 0, div_orb_e.size()) if (div_orb_e[i] < 0) ERR("Wrong with count orbit idx: " + NAV4(div_orb_e, i, ne, ns));
		for_Int(i, 0, div_orb_c.size()) if (div_orb_c[i] < 0) ERR("Wrong with count orbit idx: " + NAV4(div_orb_c, i, ne, ns));
#endif
		div_orb_es.push_back(div_orb_e);
		div_orb_cs.push_back(div_orb_c);
	}
#ifdef _ASSERTION_
	Int x(0);
	for_Int(i, 0, ns.size()) {
		for_Int(j, 0, div_orb_es[i].size()) ++x;
		for_Int(j, 0, div_orb_cs[i].size()) ++x;
	}
	if (x != SUM(ns)) ERR("Wrong with count orbit idx: " + NAV2(SUM(ns), x));
#endif
	if (div_orb_es.size() != ns.size()) ERR("Wrong with count div's orbit idx list finding: " + NAV2(div_orb_es.size(), ns.size()));
	if (div_orb_cs.size() != ns.size()) ERR("Wrong with count div's orbit idx list finding: " + NAV2(div_orb_cs.size(), ns.size())); 
	//DBG("find_div_orb() FINISHED  !!!");
	return std::make_tuple(div_orb_es, div_orb_cs);
}

Int ComDivs::find_cf_idx(const VecOnb& newcf)
{
	//DBG("find_cf_idx BEGIN :set  base" + NAV3(ne, ns, bf(ns[0], ne[0])));
	VecInt base(ns.size() + 1, 1);
	for_Int(i, 1, base.size()) base[i] = base[i - 1] * bf(ns[i - 1], ne[i - 1]);
	VecInt rep(cf.size(), 0);
	for_Int(i, 0, cf.size()) rep[i] = cf[i].idx();
	//DBG("find_cf_idx BEGIN : base, rep finished" + NAV3(base, rep, FreeBaseCode(rep, base)));2
	return FreeBaseCode(rep, base);
}

//---------------------------------------------Private function---------------------------------------------

//in this Function, the base's length = the rep's length + 1.
Int ComDivs::FreeBaseCode(VecInt& rep, VecInt& base) const
{
	Int idx(0);
	for_Int(i, 0, rep.size()) idx += base[i] * rep[i];
	if (idx >= base[base.size() - 1]) {
		MatInt test_tmp(3, ne.size(), 0);
		test_tmp[0] = ne; test_tmp[1] = ns; test_tmp[2] = rep;
		ERR(STR("read-in error with ") + NAV3(idx,test_tmp, base));
	}
	return idx;
}

// The base[0] = 1.
VecInt ComDivs::FreeBaseDeCode(Int idx, VecInt& base) const
{
	VecInt rep(base.size() - 1, 0);
	for_Int(i, 0, rep.size()) rep[i] = (idx % base[i + 1]) / base[i];
	if (idx >= base[base.size() - 1]) {
		MatInt test_tmp(3, ne.size(), 0);
		test_tmp[0] = ne; test_tmp[1] = ns; test_tmp[2] = rep;
		ERR(STR("read-in error with ") + NAV3(idx, test_tmp, base));
	}
	return rep;
}

// Only suit for fermion sign changed by c_i^+ c_j or c_j^+ c_i
Int ComDivs::count_electrons(const Int& lw, const Int& up) const
{
	Int counter(0), nele(0);

	for_Int(i, 0, ns.size()) {
		for_Int(j, 0, ns[i]) {
			if (counter >= lw && cf[i].isocc(j)) nele++;
			counter++;
			if(counter >= up)	return nele;
		}
	}
#ifdef _ASSERTION_
	if (counter > up - lw + 1) ERR("Wrong with Idxing orbit" + NAV2(lw, up));
#endif
}