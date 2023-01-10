/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022-2023
*/

#include "apizen.h"

// APIzen::APIzen(const MyMpi& mm_i, const Prmtr& prmtr_i, const Model& mdl_i):
// 	mm(mm_i), p(prmtr_i), mdl(mdl_i), iter_cnt(0), gfloc(p.norbs * 2, p),
// 	se_err({ 1.e99,1.e98,1.e97 ,1.e96 ,1.e95}), seloc(p.norbs * 2, p),
// 	//se_err({ 1.e99 }),
// 	p_n(0.),
// 	m_m(0.),
// 	d_p(0.),
// 	p_p(0.)
// {
// 	// make random seed output together
// 	{ mm.barrier(); SLEEP(1); }
// 	Bath bth(mm, p);
// 	Impurity imp(mm, p, bth, mdl);
// 	NORG norg(mm, p);

// 	norg.up_date_h0_to_solve(imp.h0);
// 	ImGreen g0imp(p.norbs * 2, p);	imp.find_g0(g0imp);						if (mm) g0imp.write("g0imp", iter_cnt);
// 	ImGreen gfimp(p.norbs * 2, p);	norg.get_g_by_KCV_spup(gfimp);			if (mm) gfimp.write("gfimp", iter_cnt);// tested.
// 	ImGreen seimp(p.norbs * 2, p);
// 	seimp = g0imp.inverse() - gfimp.inverse();								if (mm) seimp.write("seimp", iter_cnt);

// 	mdl.find_gfloc(seloc,gfloc);												if (mm) gfloc.write("gfloc");
//     // update_phy();                                                               log("se_err_converge");

//     //ReGreen g0imp_re(p.norbs * 2, p);	imp.find_g0(g0imp_re);					if (mm) g0imp_re.write("g0imp_re");
//     //ReGreen gimp_re(p.norbs * 2, p);	norg_stage2.find_g(gimp_re);			if (mm) gimp_re.write("gimp_re");
// 	//ReGreen se_re = g0imp_re.inverse() - gimp_re.inverse();					if (mm) se_re.write("se_re");
// 	//ReGreen gfloc_re(p.norbs * 2, p);	mdl.find_gfloc(se_re, gfloc_re);		if (mm) gfloc_re.write("gfloc_re");

// }





APIzen::APIzen(const MyMpi& mm_i, Prmtr& prmtr_i, const Str& file, const Int test_mode_i) :
	mm(mm_i), p(prmtr_i),
	num_orbital(prmtr_i.iqust_orbiatal), num_omg(prmtr_i.num_omg),
	imfrq_hybrid_function(num_orbital,num_omg,0.),
	solver_eimp_data(num_orbital,0.), orbitdegenerate_idx(num_orbital, 0),
	num_nondegenerate(-1), test_mode(test_mode_i)
{
	// WRN("BEGIN:: read_ZEN(file);");
	read_NORG(file);
	if (test_mode == 0) read_ZEN(file);
	// WRN("BEGIN:: fitting();");
	p.templet_restrain = restrain; p.templet_control = distribute;
	p.after_modify_prmtr();
	fitting();
	// WRN("ENDING:: fitting();");
	p.nimp = 2 * num_nondegenerate;
	p.t_ose = std::move(t_ose);
	p.t_hyb = std::move(t_hyb);
	p.mu.reset(p.nimp, 0.);
	if (test_mode) {
		for_Int(i, 0, muvec.size() * p.nimp)p.mu[i] = mune;
		p.u_hbd = Uc; p.j_ob = Jz;
	}
	else {
		for_Int(i, 0, muvec.size() * 2)p.mu[i] = muvec[int(i / 2)];
		p.u_hbd = Uc; p.j_ob = Jz;
	}
	//DBG("test for after_modify_prmtr()")
	p.after_modify_prmtr();



}

void APIzen::fitting()
{
	VEC<Bath> bthset;
	Int osesize(0), hopsize(0);
	muvec.reset(num_nondegenerate, 0.);
	for_Int(i, 0, num_nondegenerate) {
		Bath bth(mm, p);
		ImGreen hby_i(1, p);
		Int orbit_num(-1);
		for_Int(j, 0, num_orbital) {
			if (orbitdegenerate_idx[j] == i + 1) {
				orbit_num = j;
				break;
			}
		}
		for_Int(j, 0, hby_i.nomgs) hby_i.g[j] = -imfrq_hybrid_function[orbit_num][j];
		if (test_mode) { if (mm)WRN("In the test mode!!") }
		else
		{
			if (mm) hby_i.write("hby_i", orbit_num + 1);
			bth.bath_fit(hby_i, 0);
			if (mm) WRN(NAV2(bth.ose, bth.hop));
			test_for_fitting(bth, hby_i, orbit_num);
		}
		bthset.push_back(std::move(bth));
		// muvec[i] = -solver_eimp_data[orbit_num];
		muvec[i] = solver_eimp_data[orbit_num];
	}
	for_Int(i, 0, num_nondegenerate) {
		if (test_mode) {
			//VecReal tem1({ -0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4 });
			VecReal tem1(bathose);
			t_ose.push_back(tem1);
			t_ose.push_back(std::move(tem1));
			// VecReal tem2({ -4,-3,-2,-1,0,1,2,3,4 });
			VecReal tem2(bathhop);
			t_hyb.push_back(tem2);
			t_hyb.push_back(std::move(tem2));
		}
		else
		{
			t_ose.push_back(bthset[i].ose);
			t_ose.push_back(bthset[i].ose);
			t_hyb.push_back(bthset[i].hop);
			t_hyb.push_back(bthset[i].hop);
		}
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
			while (test_mode ==0){
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
