/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022-2023
*/

#include "apizen.h"

APIzen::APIzen(const MyMpi& mm_i, Prmtr& prmtr_i, const Str& file, const Int test_mode_i) :
	mm(mm_i), p(prmtr_i),num_omg(prmtr_i.num_omg),
	imfrq_hybrid_function(num_orbital,num_omg,0.),
	solver_eimp_data(num_orbital,0.), orbitdegenerate_idx(num_orbital, 0),
	num_nondegenerate(-1), test_mode(test_mode_i)
{
	read_ZEN(file);	num_orbital = norbs;
	p.templet_restrain = restrain; p.templet_control = distribute;
	p.after_modify_prmtr();
	// fitting();
	Bath bth(mm, p);
	ImGreen hb(nband, p);
	for_Int(j, 0, hb.nomgs) for_Int(i, 0, nband)  hb.g[j][i][i] = imfrq_hybrid_function[i][j];

	bth.number_bath_fit(hb, dmft_cnt, 0);
	
	p.norbs = 2 * num_nondegenerate;
	p.eimp.reset(p.norbs, 0.);
	// if (test_mode) {
	// 	for_Int(i, 0, muvec.size() * p.nimp)p.mu[i] = mune;
	// 	p.u_hbd = Uc; p.j_ob = Jz;
	// }
	// else {
	// 	for_Int(i, 0, muvec.size() * 2)p.mu[i] = muvec[int(i / 2)];
	// 	p.u_hbd = Uc; p.j_ob = Jz;
	// }
	p.hubbU = Uc;
	p.after_modify_prmtr();


	{
		// Impurity imp(mm, p, bth, mdl);
		NORG norg(mm, p);

		// norg.up_date_h0_to_solve(imp.h0);
		// ImGreen g0imp(p.norbs * 2, p);	imp.find_g0(g0imp);						if (mm) g0imp.write("g0imp", dmft_cnt);
		// ImGreen gfimp(p.norbs * 2, p);	norg.get_g_by_KCV_spup(gfimp);			if (mm) gfimp.write("gfimp", dmft_cnt);
		// ImGreen seimp(p.norbs * 2, p);	seimp=g0imp.inverse()-gfimp.inverse();	if (mm) seimp.write("seimp", dmft_cnt);
	}
}

void APIzen::test_for_fitting(const Bath& bth, const ImGreen& hby_i, Int num)
{ // test for fitting(fit_hb)
	ImGreen fit_hb(1, p);
	const VecCmplx E = cmplx(bth.ose);
	const VecCmplx V = cmplx(bth.hop);
	VecCmplx Z(p.nbath);
	for_Int(n, 0, fit_hb.nomgs)
	{
		const Cmplx z = Cmplx(0., hby_i.omg(n));
		Z = z;
		VecCmplx S = INV(E - Z);
		Cmplx hyb = SUM(V * S * V.co());
		fit_hb[n][0][0] = hyb;
	}
	fit_hb.write("hb_fitting", num + 1);
}

void APIzen::read_ZEN(const Str& file)
{
	{// hyb.in
		Str hybdata(file + ".hyb.in");
		IFS ifs(hybdata);
		if (!ifs) {
			ERR(STR("file opening failed with ") + NAV(hybdata))
		}
		else {
			Int i(0);
			bool swicher(true);
			while (true) {
				VecReal omg(num_omg, 0.), re(num_omg, 0.), im(num_omg, 0.);
				Int drop_Int(0);
				Real drop_omg(0.), drop_re(0.), drop_im(0.), drop_err1(0.), drop_err2(0.);
				if (i >= num_orbital) break;
				for_Int(j, 0, num_omg) {
					if (swicher) {
						ifs >> drop_Int;
					}
					ifs >> omg[j]; ifs >> re[j]; ifs >> im[j];
					ifs >> drop_err1; ifs >> drop_err2;
					//if(i==1)DBG(NAV2(re[j], im[j]));
					if (!ifs) ERR(STR("read_ZEN-in error with ") + NAV(hybdata));
					swicher = true;
				}
				imfrq_hybrid_function[i] = cmplx(re, im);
				while (true) {
					ifs >> drop_Int;
					if (drop_Int - 1 == i + 1) {
						swicher = false;
						break;
					}
					ifs >> drop_omg; ifs >> drop_re; ifs >> drop_im; ifs >> drop_err1; ifs >> drop_err2;
					if (!ifs) break;
				}
				++i;
			}
		}
		ifs.close();
	}

	{// eimp.in
		Str eimpdata(file + ".eimp.in");
		IFS ifs(eimpdata);
		if (!ifs) {
			ERR(STR("file opening failed with ") + NAV(eimpdata))
		}
		else {
			Int drop_Int(0);
			for_Int(i, 0, num_orbital) {
				ifs >> drop_Int;
				ifs >> solver_eimp_data[i]; ifs >> orbitdegenerate_idx[i];
				if (!ifs) ERR(STR("read_ZEN-in error with ") + NAV(eimpdata));
			}
			if (test_mode) num_nondegenerate = 1;
			else num_nondegenerate = MAX(orbitdegenerate_idx);
		}
		if (num_nondegenerate <= 0)ERR(STR("read_ZEN-in error with ") + NAV2(eimpdata, num_nondegenerate));
		ifs.close();
	}
	
	{// norg.in
		Str norgdata(file+".norg.in");
		IFS ifs(norgdata);
		if (!ifs) {
			ERR(STR("file opening failed with ") + NAV(norgdata))
		}
		else {
			Str strr;
			while (test_mode == 0){
				ifs >> strr;
				if(strr == "nband") {
					ifs >> strr;
					ifs >> nband; 
				}
				if(strr == "norbs") {
					ifs >> strr;
					ifs >> norbs; 
				}
				if(strr == "Uc") {
					ifs >> strr;
					ifs >> Uc; 
				}
				if(strr == "Jz") {
					ifs >> strr;
					ifs >> Jz; 
				}
				if(strr == "mu") {
					ifs >> strr;
					ifs >> mu; 
				}
				if (strr == "restrain") {
					Int l2, l1, r1, r2;
					ifs >> strr;
					ifs >> strr;
					ifs >> l2;	ifs >> l1; ifs >> strr; ifs >> r1; ifs >> r2;
					restrain = { 0, l2, l1, 0, r1, r2 };
					WRN("Finish the restrain input:" + NAV(restrain));
				}
				if (strr == "distribute") {
					Int l2, l1, m0, r1, r2;
					ifs >> strr;
					ifs >> strr;
					ifs >> l2;	ifs >> l1; ifs >> m0; ifs >> r1; ifs >> r2;
					distribute = { 1, l2, l1, m0, r1, r2 };
					WRN("Finish the division distribute input:" + NAV(distribute));
				}

				if (strr == "norm_mode") {
					ifs >> strr;
					ifs >> norm_mode;
					if (norm_mode)WRN("Dear cutomer,you are now in the test mode PLEASE modify the console file manually.");
				}
				if (strr == "test_mode") {
					ifs >> strr;
					ifs >> test_mode;
					if (test_mode)WRN("Dear cutomer,you are now in the test mode PLEASE modify the console file manually.");
				}
				if (strr == "fast_mode") {
					ifs >> strr;
					ifs >> fast_mode;
					if (fast_mode)ERR("Sorry, the fast mode will coming Soon...");
				}
				
				if (test_mode == 1 && strr == "bathose") {
					bathose.reset(SUM(distribute) - 1);
					Real input;
					ifs >> strr;
					for_Int(i, 0, bathose.size()) {
						ifs >> input;
						bathose[i] = input;
					}
				}
				if (test_mode == 1 && strr == "bathhop") {
					bathhop.reset(SUM(distribute) - 1);
					Real input;
					ifs >> strr;
					for_Int(i, 0, bathhop.size()) {
						ifs >> input;
						bathhop[i] = input;
					}
					WRN("Finish the test_mode input: nimp = 2" + NAV4(Uc, Jz, bathose, bathhop));
				}
				if (!ifs) break;
			}
		}
		ifs.close();
	}
}
