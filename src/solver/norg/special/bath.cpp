#include "bath.h"

Bath::Bath(const MyMpi& mm_i, const Prmtr& prmtr_i) :
	mm(mm_i), p(prmtr_i), nb(p.nI2B[0]), hb(1, p), uur(mm.id()), ose(p.nI2B[0]), hop(p.nI2B[0]), info(p.nband,5),
	vec_ose(p.nband), vec_hop(p.nband)
{
	// make random seed output together
	{ SLEEP(1); mm.barrier(); }
	// init_ose_hop();
}

void Bath::bath_fit(const ImGreen& hb_i, Int iter)
{
	read_ose_hop(); IFS ifs("ose_hop.txt");
	for_Int(band_i, 0, p.nband){
		if(p.nband != hb_i[0].nrows()) ERR("some thing wrong with the hybrid function.")
		for_Int(i, 0, hb_i.nomgs) hb[i] = hb_i[i][band_i][band_i];
		ose.reset(p.nI2B[band_i]); hop.reset(p.nI2B[band_i]); nb = p.nI2B[band_i];
		if(ifs) {ose = vec_ose[band_i]; hop = vec_hop[band_i];} 
		else init_ose_hop();
		const VecReal a0 = concat(ose, hop);
		Real err;
		VecReal a;
		Int nmin;
		std::tie(err, a, nmin) = bath_fit_contest(a0);
		ose = a.truncate(0, nb);
		hop = a.truncate(nb, nb + nb);
		regularize_ose_hop();
		if(!ifs) {vec_ose.push_back(ose);vec_hop.push_back(hop);}
		else {vec_ose[band_i] = ose; vec_hop[band_i] = hop;}
		// if(mm) WRN(NAV2(vec_ose.size(),vec_hop.size()));
		if (mm) {
			const HybErr hyberr(p, hb, nb);
			const VecReal a = concat(ose, hop);
			Real err = hyberr(a);
			Real err_crv = hyberr.err_curve(a);
			Real err_reg = hyberr.err_reg(a);
			//Real err_bsr = hyberr.err_bsr(a);
			Real a_norm = a.norm();
			using namespace std;
			cout << setw(4) << band_i+1 << "  " << NAV5(nmin, err, err_crv, err_reg,/* err_bsr,*/ a_norm) << "  " << present() << endl;
			NAV5(Int(info[band_i][0]=Real(nmin)), info[band_i][1]=err, info[band_i][2]=err_crv, info[band_i][3]=err_reg, /*err_bsr,*/ info[band_i][4]=a_norm);
		}
	}
}

void Bath::bath_fit(const ImGreen& hb_i, VecInt or_deg)
{
	read_ose_hop();IFS ifs("ose_hop.txt");
	for_Int(degi, 0, MAX(or_deg)) {
		Int count(0), orb_rep(-1);
		VecCmplx hb_fit(hb.nomgs); 
		for_Int(i, 0, hb_i.norbs) {
			count++;
			if(or_deg[i] - 1 == degi) for_Int(n, 0, hb.nomgs) hb_fit[n] += hb_i.g[n][i][i];
		}
		if(p.nband != hb_i[0].nrows()) ERR("some thing wrong with the hybrid function.")
		for_Int(i, 0, hb.nomgs) hb[i] = hb_fit[i] / count;

		for_Int(j, 0, p.norg_sets) {orb_rep = j; if(or_deg[j] == degi + 1) break;}
		ose.reset(p.nI2B[2*orb_rep]); hop.reset(p.nI2B[2*orb_rep]); nb = p.nI2B[2*orb_rep];
		if(ifs) {ose = vec_ose[2*orb_rep]; hop = vec_hop[2*orb_rep];} 
		else init_ose_hop();
		const VecReal a0 = concat(ose, hop);
		Real err;
		VecReal a;
		Int nmin;
		std::tie(err, a, nmin) = bath_fit_contest(a0);
		ose = a.truncate(0, nb);
		hop = a.truncate(nb, nb + nb);
		regularize_ose_hop();
		for_Int(i, 0, p.nband) if(or_deg[i]-1 == degi) vec_ose[i] = ose;
		for_Int(i, 0, p.nband) if(or_deg[i]-1 == degi) vec_hop[i] = hop;
		// if(mm) WRN(NAV2(vec_ose.size(),vec_hop.size()));
		if (mm) {
			const HybErr hyberr(p, hb, nb);
			const VecReal a = concat(ose, hop);
			Real err = hyberr(a);
			Real err_crv = hyberr.err_curve(a);
			Real err_reg = hyberr.err_reg(a);
			//Real err_bsr = hyberr.err_bsr(a);
			Real a_norm = a.norm();
			using namespace std;
			cout << setw(4) << degi << "  " << NAV5(nmin, err, err_crv, err_reg,/* err_bsr,*/ a_norm) << "  " << present() << endl;
			NAV5(Int(info[degi][0]=Real(nmin)), info[degi][1]=err, info[degi][2]=err_crv, info[degi][3]=err_reg, /*err_bsr,*/ info[degi][4]=a_norm);
		}
	}
}


VecReal Bath::next_initial_fitting_parameters(const VecReal& a0, const Int& ntry_fine, Int& itry)
{
	VecReal a = a0;
	++itry;
	if (itry == 1) {
	}
	else if (itry <= ntry_fine) {
		perturb(a, 0.2, 12);
	}
	else if (itry <= ntry_fine * 2) {
		perturb(a, 0.4, 8);
	}
	else if (itry <= ntry_fine * 4) {
		perturb(a, 1.0, 16);
	}
	else {
		Real exponent = 4.;
		for_Int(i, 0, a.size()) {
			Real re = uur() - 0.5;
			a[i] = SIGN(std::pow(10., (8 * ABS(re) - exponent)), a[i]);
		}
	}
	return a;
}

std::tuple<Real, VecReal, Int> Bath::bath_fit_contest(const VecReal& a0)
{
	const HybErr hyberr(p, hb, nb);
	const Int np = a0.size();
	const Int ntry_fine = MAX(16, 3 * mm.np() - 1);
	const Int ntry = MAX(128 * ntry_fine, 2000);
	const Real tol = 1.e-12;
	Int nmin = 0;		// number of fittings reaching the minimum
	MPI_Status status;
	VecReal a(np);
	VecReal a_optm = a0;
	Real err_optm = hyberr(a_optm);

	if (mm) {
		Int ntot = ntry;
//		Int ntot = 2;
		Int nsnd = 0;
		Int nrcv = 0;
		Int nfst = 1;
		Int itry = 0;
		while (nrcv < ntot) {
			if (nfst < mm.np()) {
				if (nsnd < ntot) {
					a = next_initial_fitting_parameters(a0, ntry_fine, itry);
					mm.Send(a, nfst, 1);
					++nsnd;
				}
				else {
					mm.Send(a, nfst, 0);
				}
				++nfst;
			}
			else {
				mm.Recv(a, status);
				++nrcv;
				Int sndr = status.MPI_SOURCE;
				Real err = hyberr(a);
				if (false) {
					Real err_crv = hyberr.err_curve(a);
					Real err_reg = hyberr.err_reg(a);
					//Real err_bsr = hyberr.err_bsr(a);
					Real a_norm = a.norm();
					WRN(NAV5(sndr, ntot, nsnd, nrcv, itry) + ", " + NAV5(err_optm, err, err_crv, err_reg,/* err_bsr,*/ a_norm));
				}
				if (err_optm - tol > err) { nmin = 0; }
				if (err_optm > err) { a_optm = a; err_optm = err; }
				nmin += err - err_optm < tol;

				if (nsnd < ntot) {
					a = next_initial_fitting_parameters(a0, ntry_fine, itry);
					mm.Send(a, sndr, 1);
					++nsnd;
				}
				else {
					mm.Send(a, sndr, 0);
				}
			}
		}
	}
	else {
		while (true) {
			mm.Recv(a, status, mm.ms());
			if (status.MPI_TAG == 0) break;
			FitMrq<HybErr> mrq(hyberr.x, hyberr.y, hyberr.sig, a, hyberr, tol);
			Int mrq_fit_info = mrq.fit();
			mm.Send(mrq.a, mm.ms(), 1);
		}
	}

	mm.Bcast(err_optm);
	mm.Bcast(a_optm);
	mm.Bcast(nmin);
	return std::make_tuple(err_optm, a_optm, nmin);
}
//rely on imp model's frame
MatReal Bath::find_hop() const
{
    MatReal h0(0, 0, 0.);
    for_Int(i, 0, p.nband) {
        const Int nb_i = p.nI2B[i * 2];
        MatReal h0_i(1 + nb_i, 1 + nb_i, 0.);
        for_Int(j, 0, nb_i) {
            h0_i[0][j + 1] = vec_hop[i][j];
            h0_i[j + 1][0] = vec_hop[i][j];
            h0_i[j + 1][j + 1] = vec_ose[i][j];
        }
        h0_i.reset(direct_sum(h0_i, h0_i));
        h0.reset(direct_sum(h0, h0_i));
    }
    return h0;
}

void Bath::write_ose_hop(Int iter_cnt) const {
	using namespace std;
	// OFS ofs_app_ose(p.ofx + ".ose.txt", std::ios::app);
	// OFS ofs_ose;ofs_ose.open("ose.txt");
	OFS ofs;ofs.open("ose_hop.txt");
	for_Int(band_i, 0, p.nband)	{
		ofs << iofmt("sci");
		ofs << setw(4) << iter_cnt << setw(4) << band_i;
		for_Int(i, 0, p.nI2B[2 * band_i]) { ofs << "\t" << setw(w_Real) << vec_ose[band_i][i]; }
		ofs << endl;
	}
	for_Int(band_i, 0, p.nband)	{
		ofs << iofmt("sci");
		ofs << setw(4) << iter_cnt << setw(4) << band_i;
		for_Int(i, 0, p.nI2B[2 * band_i]) { ofs << "\t" << setw(w_Real) << vec_hop[band_i][i]; }
		ofs << endl;
	}
	ofs << endl;ofs << endl;
}

void Bath::read_ose_hop(Int iter_cnt) {
	using namespace std;
	// IFS ifs_ose;ifs_ose.open("ose.txt");
	IFS ifs("ose_hop.txt");
	Str strr;
	if(ifs)for_Int(band_i, 0, p.nband)	{
		ifs >> strr; ifs >> strr;
		VecReal ose_t(p.nI2B[2 * band_i], 0);
		for_Int(i, 0, p.nI2B[2 * band_i]) { ifs >> ose_t[i]; }
		// vec_ose.push_back(ose_t);
		vec_ose[band_i] = ose_t;
		// if(mm) WRN(NAV(ose_t));
	}		
	if(ifs)for_Int(band_i, 0, p.nband)	{
		ifs >> strr; ifs >> strr;
		VecReal hop_t(p.nI2B[2 * band_i], 0);
		for_Int(i, 0, p.nI2B[2 * band_i]) { ifs >> hop_t[i]; }
		// vec_hop.push_back(hop_t);
		vec_hop[band_i] = hop_t;
		// if(mm) WRN(NAV(hop_t));
	}
}

//------------------------------------------------------------------ print out ------------------------------------------------------------------
void Bath::regularize_ose_hop() {
	slctsort(ose, hop);
	// after a unitary transformation, hop can always be
	// non-negative when there is only one impurity site
	hop = ABS(hop);
	for_Int(i, 0, ose.size()) if (ABS(ose[i]) > p.hubbU * 8) hop[i] = ose[i] = 0.;
}