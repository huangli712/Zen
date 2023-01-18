#include "bath.h"

Bath::Bath(const MyMpi& mm_i, const Prmtr& prmtr_i) :
	mm(mm_i), p(prmtr_i), nb(p.nI2B[0]), hb(1, p), uur(mm.id()), ose(p.nI2B[0]), hop(p.nI2B[0])
{
	// make random seed output together
	{ SLEEP(1); mm.barrier(); }
	// init_ose_hop();
}

void Bath::bath_fit(const ImGreen& hb_i, Int iter)
{
	for_Int(band_i, 0, p.nband){
		if(p.nband != hb_i[0].nrows()) ERR("some thing wrong with the hybrid function.")
		for_Int(i, 0, hb_i.nomgs) hb[i] = hb_i[i][band_i][band_i];
		ose.reset(p.nI2B[band_i]); hop.reset(p.nI2B[band_i]); nb = p.nI2B[band_i];
		init_ose_hop();
		const VecReal a0 = concat(ose, hop);
		Real err;
		VecReal a;
		Int nmin;
		std::tie(err, a, nmin) = bath_fit_contest(a0);
		ose = a.truncate(0, nb);
		hop = a.truncate(nb, nb + nb);
		regularize_ose_hop();
		vec_ose.push_back(ose);vec_hop.push_back(hop);
		if (mm) {
			const HybErr hyberr(p, hb, nb);
			const VecReal a = concat(ose, hop);
			Real err = hyberr(a);
			Real err_crv = hyberr.err_curve(a);
			Real err_reg = hyberr.err_reg(a);
			//Real err_bsr = hyberr.err_bsr(a);
			Real a_norm = a.norm();
			using namespace std;
			cout << setw(4) << iter << "  " << NAV6(band_i,nmin, err, err_crv, err_reg, /*err_bsr,*/ a_norm) << "  " << present() << endl;
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
	const Int ntry_fine = MAX(16, mm.np());
	const Int ntry = MAX(8 * ntry_fine, 2000);
	const Real tol = 1.e-10;
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